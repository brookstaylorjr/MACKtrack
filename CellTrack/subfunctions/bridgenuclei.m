function [label_out] = bridgenuclei(subobj_in,obj_cc, cutoff,verbose)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [label_out] = bridgenuclei(label_in,cutoff,verbose)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% BRIDGENUCLEI checks all objects in a label mask, seeing if any can be bridged to yield a 
% morphologically-defined nucleus.
% [Note: objects need to be touching each other in order to be bridged.]
%
% subobj_in     label matrix ('subobjects') to be merged
% obj_cc        bwconncomp structure of completely grouped objects
% cutoffs       morphological cutoffs: determine whether object is kept/merged/dropped. Cutoffs
%               must contain: area ([min max]), eccentricty (min), and compactness ([low high])
% verbose       boolean, 1 displays merging/testing output
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%%
% Input checking - return if empty matrix is passed
if nargin <4
    verbose=0;
end
if length(unique(subobj_in))<2
    label_out=double(subobj_in);
    return;
end

% Dilate over of watershed borders, compute subobject and grouped-object CC structure
cc_in = label2cc(subobj_in,0);
obj_rprops = cell2mat(struct2cell(regionprops(obj_cc,'Area','Solidity','Perimeter')))';
obj_rprops = [obj_rprops(:,1),obj_rprops(:,3),obj_rprops(:,2)]; % (Rearrange to put solidity @ end)
obj_rprops = mat2cell(obj_rprops,ones(1,size(obj_rprops,1)),3);

% Assign groups of subobjects to corresponding parent object (in obj_cc)
get_obj = @(pxlist) unique(subobj_in(pxlist));
obj_match = cellfun(get_obj, obj_cc.PixelIdxList,'UniformOutput',0);
obj_groups = cell(size(cc_in.PixelIdxList));
for i = 1:length(obj_match)
    obj_groups(obj_match{i}) = obj_match(i);
end

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
                num2str((obj_rprops{i}(2).^2)/(4*pi*obj_rprops{i}(1))),',  s = ',num2str(obj_rprops{i}(3)),')'];
        if pass_1(i)
            disp(['I. ADDED [',num2str(sum(pass_1(1:i))),']: ', str_1]);
        else
            disp(['FAILED ', str_1]);
        end
    end
end
remaining_obj = find(~pass_1); % Remaining grouped objects

if ~isequal(size(obj_match),size(obj_rprops))
    obj_match = obj_match';
end
% PASS 2: Soft pass (with solidity check) on (combined) objects with 1-2 subobjects
pass_2 = (cellfun(@length,obj_match(remaining_obj))<=2) & cellfun(soft_plus,obj_rprops(remaining_obj));
if verbose
    add_list = remaining_obj(pass_2);
    for i = 1:length(add_list)
        disp(['III. ADDED [',num2str(i+length(pixelidx_out)),']: obj containing [',num2str(obj_match{add_list(i)}'),']. (a = ',...
            num2str(obj_rprops{add_list(i)}(1)),', c = ',...
            num2str((obj_rprops{add_list(i)}(2).^2)/(4*pi*obj_rprops{add_list(i)}(1))),')'])
    end
end
pixelidx_out = [pixelidx_out(:); reshape(obj_cc.PixelIdxList(remaining_obj(pass_2)),sum(pass_2),1)];
remaining_obj = remaining_obj(~pass_2);

% Drop out remaining single-subobject shapes to prevent further analysis
remaining_obj(cellfun(@length,obj_match(remaining_obj))==1) = [];

% PASS 3: Strong pass, trying to combine subsets of subobjects
remaining_subobj = cell2mat(obj_match(remaining_obj));
cc1 = cc_in;
counter = 0;
while ~isempty(remaining_subobj) && (counter<100)
    [cc1, remaining_subobj] = mergeNeighbors(cc1, remaining_subobj, obj_groups, cutoff, 'hard',3, verbose);
    counter = counter+1;
end
remaining_subobj = cell2mat(obj_match(remaining_obj)); % Get original list of subobjects (minus newly combined ones)  
remaining_subobj(ismember(remaining_subobj,find(cellfun(@isempty,cc1.PixelIdxList)))) = [];

% PASS 4: Of remaining objects, see if any is composed of subobjects that all hard_pass
label1 = labelmatrix(cc1);
label1(~ismember(label1,remaining_subobj)) = 0;
obj_cc1 = bwconncomp(label1>0,4);
if obj_cc1.NumObjects>0
    rprops1 = cell2mat(struct2cell(regionprops(label1,'Area','Perimeter')))';
    rprops1 = mat2cell(rprops1,ones(1,size(rprops1,1)),2);
    get_obj = @(pxlist) unique(label1(pxlist));
    obj_match = cellfun(get_obj, obj_cc1.PixelIdxList,'UniformOutput',0)';
    pass_orig = cellfun(hard_pass, rprops1);
    all_pass = @(subobj) min(pass_orig(subobj));
    pass_4 = cellfun(all_pass, obj_match);
    if verbose
        add_list = obj_match(pass_4);
        for i = 1:length(add_list)
            disp(['II. ADDED [',num2str(i+length(pixelidx_out)),']: obj containing (all strict-passing) subobj [',...
                num2str(add_list{i}'),']'])
        end
    end
    pixelidx_out = [pixelidx_out(:); cc1.PixelIdxList(cell2mat(obj_match(pass_4)))];
    remaining_subobj(ismember(remaining_subobj,cell2mat(obj_match(pass_4)))) = [];
end

% PASS 5: Try to combine subobjects (up to 4) and pass with "soft" criteria
if length(remaining_subobj)>0
    [cc1, remaining_subobj] = mergeNeighbors(cc1, remaining_subobj, obj_groups, cutoff,'soft',4, verbose);
end

% PASS 6: Perform final soft-pass on remaining (uncombined) objects
label1 = labelmatrix(cc1);
label1(cell2mat(pixelidx_out)) = 0;
cc2 = label2cc(label1);
if cc2.NumObjects>0
    rprops_new = cell2mat(struct2cell(regionprops(cc2,'Area','Perimeter')))';
    rprops_new = mat2cell(rprops_new,ones(1,size(rprops_new,1)),2);
    pass_5 = find(cellfun(soft_pass, rprops_new));
else
    pass_5 = [];
end

% Make output label matrix
cc_out.PixelIdxList = [pixelidx_out;cc2.PixelIdxList(pass_5)];
cc_out.NumObjects = length(cc_out.PixelIdxList);
cc_out.ImageSize = size(subobj_in);
cc_out.Connectivity = 4;


% Make label_out
label_out = double(labelmatrix(cc_out));



function [cc_out, subobj_out] = mergeNeighbors(cc_in, subobj, obj_groups, cutoff, pass_type, dim, verbose)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MERGENEIGHBORS finds connected object subsets
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SETUP: calculate base regionprops
label_in = labelmatrix(cc_in);
rprops_in = cell2mat(struct2cell(regionprops(label_in,'Area','Solidity','Perimeter')))';
rprops_in = [rprops_in(:,1),rprops_in(:,3),rprops_in(:,2)]; % (Rearrange to put solidity @ end)
rprops_in = mat2cell(rprops_in,ones(1,size(rprops_in,1)),3);
cc_out = cc_in;

% Calculate neighbors
label_dil = imdilate(label_in,diskstrel(1));
tmp = label_in; tmp(tmp==0) = max(tmp(:))+1;
label_erd = imerode(tmp,diskstrel(1));
label_erd2 = imerode(label_in,diskstrel(1));
get_neighbors = @(set) unique([label_dil(set); label_erd(set); label_erd2(set)]);
neighbors_in = cellfun(get_neighbors,cc_in.PixelIdxList,'UniformOutput',0);
filter_neighbors = @(set1, set2) set1(ismember(set1,set2));
neighbors_in = cellfun(filter_neighbors,neighbors_in,obj_groups,'UniformOutput',0);

% Create morphological tests, evaluate existing objects
hard_pass = @(x) (x(1) < cutoff.Area(2)) && (x(1) > cutoff.Area(1)) && (((x(2).^2)/(4*pi*x(1)))<cutoff.Compactness(1))...
    && (x(3)>cutoff.Solidity);
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


rprops1 = cell2mat(struct2cell(regionprops(tmp_cc,'Area','Solidity','Perimeter')))';
rprops1 = [rprops1(:,1),rprops1(:,3),rprops1(:,2)]; % (Rearrange to put solidity @ end
rprops1 = mat2cell(rprops1,ones(1,size(rprops1,1)),3);
switch pass_type
    case 'hard'
        pass_fcn = hard_pass;
    case 'soft'
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
            disp(['V. MERGED subobj [',num2str(combo),']. (',tmp_str,')'])
        end      
    end
end



function [subset_list] = getSubsets(obj, all_neighbors, n, dropObj)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% GETSUBSETS finds connected subobject subsets (e.g. cobminations of touching objects)
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



