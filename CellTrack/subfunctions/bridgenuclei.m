function [label_out] = bridgenuclei(label_in,cutoff,verbose)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [label_out] = bridgenuclei(label_in,cutoff,verbose)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% BRIDGENUCLEI checks all objects in a label mask, seeing if any can be bridged to yield a 
% morphologically-defined nucleus.
% [Note: objects need to be touching each other in order to be bridged.]
%
% old_label     label matrix to be merged
% cutoffs       morphological cutoffs: determine whether object is kept/merged/dropped. Cutoffs
%               must contain: area ([min max]), eccentricty (min), and compactness ([low high])
% verbose       boolean, 1 displays merging/testing output
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if nargin <3
    verbose=0;
end

% Dilate label matrix over watershed borders - if there are no objects to bridge, then break out of function.
label_in = imclose(label_in,ones(2));
if length(unique(label_in))<2
    label_out=label_in;
    return;
end

% Calculate regionprops for base objects
in_rprops = cell2mat(struct2cell(regionprops(label_in,'Area','Perimeter')))';
in_rprops = mat2cell(in_rprops,ones(1,size(in_rprops,1)),2);

% Combine all connected objects and calculate regionprops
obj_cc = bwconncomp(label_in>0,4);
obj_rprops = cell2mat(struct2cell(regionprops(obj_cc,'Area','Solidity','Perimeter')))';
obj_rprops = [obj_rprops(:,1),obj_rprops(:,3),obj_rprops(:,2)]; % (Rearrange to put solidity @ end
obj_rprops = mat2cell(obj_rprops,ones(1,size(obj_rprops,1)),3);

% Assign initial objects (in label in) to groups of objects (in obj_cc)
get_obj = @(pxlist) unique(label_in(pxlist));
obj_match = cellfun(get_obj, obj_cc.PixelIdxList,'UniformOutput',0)';

% Create base checking functions
hard_pass = @(x) (x(1) < cutoff.Area(2)) && (x(1) > cutoff.Area(1)) && (((x(2).^2)/(4*pi*x(1)))<cutoff.Compactness(1));
soft_pass = @(x) (x(1) < cutoff.Area(2)) && (x(1) > cutoff.Area(1)) && (((x(2).^2)/(4*pi*x(1)))<cutoff.Compactness(2));
hard_plus = @(x) (x(1) < cutoff.Area(2)) && (x(1) > cutoff.Area(1)) && (((x(2).^2)/(4*pi*x(1)))<cutoff.Compactness(1))...
    && (x(3)>cutoff.Solidity);

soft_plus = @(x) (x(1) < cutoff.Area(2)) && (x(1) > cutoff.Area(1)) && (((x(2).^2)/(4*pi*x(1)))<cutoff.Compactness(2))...
    && (x(3)>cutoff.Solidity);

% PASS 1: Check intially-combined objects -> add in to pixels_out, drop from obj_cc
pass_1 = cellfun(hard_plus,obj_rprops);
pixelidx_out = obj_cc.PixelIdxList(pass_1);
if verbose
    for i = 1:length(pass_1)
        str_1 = ['obj containing [',num2str(obj_match{i}'),']. (a = ',num2str(obj_rprops{i}(1)),', c = ',...
                num2str((obj_rprops{i}(2).^2)/(4*pi*obj_rprops{i}(1))),')'];
        if pass_1(i)
            disp(['I. ADDED [',num2str(sum(pass_1(1:i))),']: ', str_1]);
        else
            disp(['FAILED ', str_1]);
        end
    end
end
remaining_obj = find(~pass_1);


% PASS 2: Of remaining (combined) objects, see if any is composed of subobjects that all hard_pass
pass_orig = cellfun(hard_pass, in_rprops);
all_pass = @(subobj) min(pass_orig(subobj));
pass_2 = cellfun(all_pass, obj_match(remaining_obj));
cc_in = label2cc(label_in,0);
if verbose
    add_list = cell2mat(obj_match(remaining_obj(pass_2)));
    for i = 1:length(add_list)
        disp(['II. ADDED [',num2str(i+length(pixelidx_out)),']: obj containing (all strict-passing) subobj [',...
            num2str(add_list(i)'),']'])
    end
end
pixelidx_out = [pixelidx_out(:); cc_in.PixelIdxList(cell2mat(obj_match(remaining_obj(pass_2))))];
remaining_obj = remaining_obj(~pass_2);


% PASS 3: soft (with solidity) pass on single/double objects
pass_3 = (cellfun(@length,obj_match(remaining_obj))<=2) & cellfun(soft_plus,obj_rprops(remaining_obj));
if verbose
    add_list = remaining_obj(pass_3);
    for i = 1:length(add_list)
        disp(['III. ADDED [',num2str(i+length(pixelidx_out)),']: obj containing [',num2str(obj_match{add_list(i)}'),']. (a = ',...
            num2str(obj_rprops{add_list(i)}(1)),', c = ',...
            num2str((obj_rprops{add_list(i)}(2).^2)/(4*pi*obj_rprops{add_list(i)}(1))),')'])
    end
end
pixelidx_out = [pixelidx_out; reshape(obj_cc.PixelIdxList(remaining_obj(pass_3)),sum(pass_3),1)];
remaining_obj = remaining_obj(~pass_3);

% Drop out remaining single/double-subobject shapes
remaining_obj(cellfun(@length,obj_match(remaining_obj))==1) = [];


% PASS 4A: Strong pass, sets of 2 or 3 subobjects
subobj = cell2mat(obj_match(remaining_obj));
while length(subobj)>0 
    [cc_in, subobj] = mergeNeighbors(cc_in, subobj, cutoff,'hard',3, verbose);
end

% PASS 4B: Two combination attempts with soft pass criteria (+solidity)
subobj = cell2mat(obj_match(remaining_obj)); % Get original list of subobjects ( minus newly combined ones)  
subobj(ismember(subobj,find(cellfun(@isempty,cc_in.PixelIdxList)))) = [];
if length(subobj)>0
    [cc_in, subobj] = mergeNeighbors(cc_in, subobj, cutoff,'soft',4, verbose);
end

% Get remaining uncombined objects that also soft-pass, make label_out
rprops_new = cell2mat(struct2cell(regionprops(labelmatrix(cc_in),'Area','Perimeter')))';
rprops_new = mat2cell(rprops_new,ones(1,size(rprops_new,1)),2);
keepers = find(cellfun(soft_pass, rprops_new));
keepers(~ismember(keepers,cell2mat(obj_match(remaining_obj)))) = [];
cc_out.PixelIdxList = [pixelidx_out;cc_in.PixelIdxList(keepers)];
cc_out.NumObjects = length(cc_out.PixelIdxList);
cc_out.ImageSize = size(label_in);
cc_out.Connectivity = 4;
label_out = double(labelmatrix(cc_out));



function [cc_out, subobj_out] = mergeNeighbors(cc_in, subobj, cutoff,pass_type, dim, verbose)


% SETUP: calculate base regionprops
label_in = labelmatrix(cc_in);
rprops_in = cell2mat(struct2cell(regionprops(label_in,'Area','Perimeter')))';
rprops_in = mat2cell(rprops_in,ones(1,size(rprops_in,1)),2);
cc_out = cc_in;

% Calculate neighbors
label_dil = imdilate(label_in,diskstrel(1));
tmp = label_in; tmp(tmp==0) = max(tmp(:))+1;
label_erd = imerode(tmp,diskstrel(1));
label_erd2 = imerode(label_in,diskstrel(1));
get_neighbors = @(set) unique([label_dil(set); label_erd(set); label_erd2(set)]);
neighbors_in = cellfun(get_neighbors,cc_in.PixelIdxList,'UniformOutput',0);


% Create morphological tests, evaluate existing objects
hard_pass = @(x) (x(1) < cutoff.Area(2)) && (x(1) > cutoff.Area(1)) && (((x(2).^2)/(4*pi*x(1)))<cutoff.Compactness(1));
soft_pass = @(x) (x(1) < cutoff.Area(2)) && (x(1) > cutoff.Area(1)) && (((x(2).^2)/(4*pi*x(1)))<cutoff.Compactness(2));
soft_plus = @(x) (x(1) < cutoff.Area(2)) && (x(1) > cutoff.Area(1)) && (((x(2).^2)/(4*pi*x(1)))<cutoff.Compactness(2))...
    && (x(3)>cutoff.Solidity);

only_hard =  cellfun(hard_pass, rprops_in);
only_fail = ~cellfun(soft_pass, rprops_in);
only_soft = ~only_fail & ~only_hard;


% 1) SEARCH all UNIQUE groupings of size <=n subobjects
groups = cell(dim-1,1);
for i = 1:length(subobj)
    tmp = getSubsets(double(subobj(i)),neighbors_in,dim);
    for j = 1:length(groups)
        groups{j} = cat(1,groups{j},tmp{j});
    end
end
tmp_cc.Connectivity = 4;
tmp_cc.ImageSize = cc_in.ImageSize;
tmp_cc.PixelIdxList = {};
tmp_cc.NumObjects = 0;
for j = 1:length(groups)
    groups{j} = unique(sort(groups{j},2),'rows');
    tmp_cc.NumObjects = tmp_cc.NumObjects+size(groups{j},1);
end

% If there's nothing to merge, break out
if tmp_cc.NumObjects==0
    cc_out = cc_in;
    subobj_out = [];
    return
end


% MEASURE regionprops for each of these possible groupings, calculate corresponding hard/soft pass info
for i = 1:length(groups)
    for j = 1:size(groups{i},1)
        tmp_cc.PixelIdxList = cat(1, tmp_cc.PixelIdxList,cell2mat(cc_in.PixelIdxList(groups{i}(j,:))));
    end
end

 
switch pass_type
    case 'hard'
        rprops1 = cell2mat(struct2cell(regionprops(tmp_cc,'Area','Perimeter')))';
        rprops1 = mat2cell(rprops1,ones(1,size(rprops1,1)),2);
        pass_fcn = hard_pass;
    case 'soft'
        rprops1 = cell2mat(struct2cell(regionprops(tmp_cc,'Area','Solidity','Perimeter')))';
        rprops1 = [rprops1(:,1),rprops1(:,3),rprops1(:,2)]; % (Rearrange to put solidity @ end
        rprops1 = mat2cell(rprops1,ones(1,size(rprops1,1)),3);
        pass_fcn = soft_plus;
end
main_pass = cellfun(pass_fcn,rprops1);
    

% Calculate the objective function savings achieved by each combination
savings = [];
map1 = 0;
for i = 1:length(groups)
    map1 = cat(1,map1,size(groups{i},1));
    if size(groups{i},1)>1
        savings = cat(1,savings,sum(only_hard(groups{i})+2*only_soft(groups{i})+4*only_fail(groups{i}),2));
    elseif size(groups{i},1)>0
        savings = cat(1,savings,sum(only_hard(groups{i})+2*only_soft(groups{i})+4*only_fail(groups{i})));
    end
end
map1 = cumsum(map1);


% RESOLVE subobject groupings in order of their savings (hard passes only, for now)
[~,order]  =sort(savings,'descend');
tmp1 = main_pass(order);
group_rows = order(tmp1);

if isempty(group_rows)
    cc_out = cc_in;
    subobj_out = [];
    return;
end
get_combo = @(idx) groups{find(map1<idx,1,'last')}(idx-map1(find(map1<idx,1,'last')),:);

if verbose
    props = cell2mat(rprops1);
    props = props(group_rows,:);
end

% Resolution: 
% a) update pixelidxlist in new_cc (combine pixels of subobjects)
% b) update neighbors
% c) record implicated subobjects -> each can only be used once.
combined_objects = [];
subobj_out = []; % Restrict subobjects to newly created clusters, since others are evaluated on first pass

for i = 1:length(group_rows)
    combo = get_combo(group_rows(i));
    if max(ismember(combo,combined_objects))==0
        combined_objects = cat(2,combined_objects,combo);
        subobj_out = cat(2,subobj_out,combo(1));
        cc_out.PixelIdxList{combo(1)} = cell2mat(cc_out.PixelIdxList(combo));
        cc_out.PixelIdxList(combo(2:end)) = cell(length(combo)-1, 1);
        if verbose
            tmp_num = props(i,:);
            tmp_str = ['a = ', num2str(tmp_num(1)),', c = ', num2str((tmp_num(2).^2)/(4*pi*tmp_num(1)))];
            if length(tmp_num)>2
                tmp_str = [tmp_str,', s = ',num2str(tmp_num(3))];
            end
            disp(['IV. COMBINED subobj [',num2str(combo),']. (',tmp_str,')'])
        end      
    end
end



function [subset_list] = getSubsets(obj, all_neighbors, n, dropObj)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% GETSUBSETS finds connected object subsets
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if nargin<4
    dropObj = [];
end

neighb1 = filterNeighbors(all_neighbors{obj},obj,dropObj);

% Get all sets centered on object - store in arrays (i.e. by size)
subset_list = cell(n-1,1);
for i = (n-1):-1:1
    if length(neighb1) >= i
        new_subsets = nchoosek(neighb1,i);
        new_subsets =  sort([obj*ones(size(new_subsets,1),1), new_subsets],2);
        subset_list{i} = new_subsets;
    end
end

% Combine subset lists into a master cell array
check_subset = {};
for i = 1:length(subset_list)
    check_subset = cat(1,check_subset, num2cell(subset_list{i},2));
end


% Work down list of connected objects - merge lists of nieghbors, look for new possible subsets from merged list
while ~isempty(check_subset)
    if length(check_subset{1}) < n
        merged_neighbors = cell2mat(all_neighbors(check_subset{1}));
        merged_neighbors = filterNeighbors(merged_neighbors,check_subset{1},dropObj);
        % Find valid subsets using (merged) neighbors of current cluster
        for j = (n-i):-1:1
            if length(merged_neighbors) >= j
                new_subsets = nchoosek(merged_neighbors,j);
                new_subsets =  sort([repmat(check_subset{1},size(new_subsets,1),1), new_subsets],2);
                curr_idx = size(new_subsets,2) - 1;
                if isempty(subset_list{curr_idx})
                    subset_list{curr_idx} = new_subsets;
                else
                    subset_list{curr_idx} = cat(1,subset_list{curr_idx},new_subsets);
                end
                
                % Drop non-unique rows out of the subset list
                [~,unq_idx] = unique(subset_list{curr_idx},'rows');
                subset_list{curr_idx} = subset_list{curr_idx}(unq_idx,:);
                % Add new objects into running list of subsets to check further
                check_subset = cat(1,check_subset,num2cell(subset_list{curr_idx},2));
            end
        end
        check_subset(1) = [];

    else % Invalid subset - drop from list
        check_subset(1)= [];
    end
end
function filtered_set = filterNeighbors(set,self,dropObj)
set(ismember(set,self)) = [];
set(set==0) = [];
set(ismember(set,dropObj)) = [];
filtered_set = unique(set);



