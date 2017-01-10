function [label_out] = bridgenuclei(subobj_in,obj_cc, cutoff, shapedef, verbose)
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
% shapedef      'Compactness' or 'Solidity' (criteria used to assess 
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
%%
% Define regionprops/filtering functions, get regionprops of objects
fcn = definefilters(cutoff,shapedef);
cc_in = label2cc(subobj_in,0);
obj_rprops = fcn.get_rprops(obj_cc);


% Assign groups of subobjects to corresponding parent object (in obj_cc)
get_obj = @(pxlist) unique(subobj_in(pxlist));
obj_match = cellfun(get_obj, obj_cc.PixelIdxList,'UniformOutput',0);
rm_zeros = @(pxlist) pxlist(pxlist>0); % Filter "improper objects" out of there.
obj_match = cellfun(rm_zeros, obj_match,'UniformOutput',0);
obj_match = obj_match(:);
obj_groups = cell(size(cc_in.PixelIdxList));
for i = 1:length(obj_match)
    obj_groups(obj_match{i}) = obj_match(i);
end


% PASS 1: Check intially-combined objects -> add in to pixels_out, drop from obj_cc
pass_1 = fcn.soft_pass(obj_rprops);
pixelidx_out = obj_cc.PixelIdxList(pass_1);
pixelidx_out = pixelidx_out(:);
if verbose
    for i = 1:length(pass_1)
        str_1 = [' containing subobj. [ ',num2str(obj_match{i}'),' ]. (Area: ',num2str(obj_rprops(i,1)),...
            ', ',shapedef, ': ',num2str(fcn.shape(obj_rprops(i,:))),')'];
        if pass_1(i)
            disp(['I. ADDED # ',num2str(sum(pass_1(1:i))),': ', str_1]);
        else
            disp(['I. FAILED obj ', str_1]);
        end
    end
    disp('- - - - - - - - - - - - - - - - - - ')
end
remaining_obj = find(~pass_1); % Remaining grouped objects



%% PASS 2: Strong pass, trying to combine subsets of subobjects
remaining_subobj = cell2mat(obj_match(remaining_obj));
cc1 = cc_in;
counter = 0;
while ~isempty(remaining_subobj) && (counter<100)
    [cc1, remaining_subobj] = mergeNeighbors(cc1, remaining_subobj, obj_groups, cutoff, shapedef, 'hard', 3, verbose);
    counter = counter+1;
end
remaining_subobj = cell2mat(obj_match(remaining_obj)); % Get original list of subobjects (minus newly combined ones)  
remaining_subobj(ismember(remaining_subobj,find(cellfun(@isempty,cc1.PixelIdxList)))) = [];


%% PASS III: Perform final area-only pass on remaining (uncombined) objects
label1 = labelmatrix(cc1);
label1(cell2mat(pixelidx_out)) = 0;
cc2 = label2cc(label1,0);
if cc2.NumObjects>0
    rprops_new = fcn.get_rprops(label1);
    pass_5 = find(fcn.end_pass(rprops_new)>0);
    if verbose
        for i = 1:length(pass_5)
            str_1 = ['(Area: ',num2str(rprops_new(pass_5(i),1)),...
                ', ',shapedef, ': ',num2str(fcn.shape(rprops_new(pass_5(i),:))),')'];
                disp(['III. ADDED # ',num2str(length(pixelidx_out)+i),': ', str_1]);
        end
    end

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


%%
function [cc_out, subobj_out] = mergeNeighbors(cc_in, subobj, obj_groups, cutoff, shapedef, pass_type, dim, verbose)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MERGENEIGHBORS finds connected object subsets
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SETUP: calculate base regionprops, get definition p
fcn = definefilters(cutoff,shapedef);
label_in = labelmatrix(cc_in);
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

% If there's nothing to merge, break out of fcn
if tmp_cc.NumObjects==0
    cc_out = cc_in;
    subobj_out = [];
    return
end

% Do morphological tests on base objects
rprops_in = fcn.get_rprops(labelmatrix(cc_in));
only_hard =  fcn.hard_pass(rprops_in);
only_fail = ~fcn.soft_pass(rprops_in);
only_soft = ~only_fail & ~only_hard;


% MEASURE regionprops for each of these possible groupings, calculate corresponding hard/soft pass info
for i = 1:length(groups)
    for j = 1:size(groups{i},1)
        tmp_cc.PixelIdxList = cat(1, tmp_cc.PixelIdxList,cell2mat(cc_in.PixelIdxList(groups{i}(j,:))));
    end
end

rprops1 = fcn.get_rprops(tmp_cc);
switch pass_type
    case 'hard'
        pass_fcn = fcn.hard_pass;
    case 'soft'
        pass_fcn = fcn.soft_pass;
end
main_pass = pass_fcn(rprops1);
    

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
    props = rprops1;
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
            tmp_str = ['Area: ', num2str(tmp_num(1)),', ',shapedef, ': ',num2str(fcn.shape(tmp_num))];
            disp(['II. MERGED subobj [ ',num2str(combo),' ]. (',tmp_str,')'])
        end      
    end
end


%%
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

%%
function filtered_set = filterNeighbors(set,self,dropObj)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% FILTERNEIGHBORS cleans up a potential set of neighboring objects
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
set(ismember(set,self)) = [];
set(set<=0) = [];
set(ismember(set,dropObj)) = [];
filtered_set = unique(set);



%%
function fcn = definefilters(cutoff,shapedef)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% DEFINEFILTERS defines anonymous functions used to get regionproperties and filter nuclear shapes
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if strcmpi(shapedef,'Compactness')
   fcn.get_rprops = @(cc) cell2mat(struct2cell(regionprops(cc,'Area','Perimeter')))';
   fcn.shape = @(x) (x(:,2).^2)./(4*pi*x(:,1));
   fcn.hard_pass = @(x) (x(:,1) < cutoff.Area(2)) & (x(:,1) > cutoff.Area(1)) & (fcn.shape(x)<cutoff.Compactness(1));
   fcn.soft_pass = @(x) (x(:,1) < cutoff.Area(2)) & (x(:,1) > cutoff.Area(1)) & (fcn.shape(x)<cutoff.Compactness(2));
else % (Solidity)
   fcn.shape = @(x) x(:,2);
   fcn.get_rprops = @(cc) cell2mat(struct2cell(regionprops(cc,'Area','Solidity')))';
   fcn.hard_pass = @(x) (x(:,1) < cutoff.Area(2)) & (x(:,1) > cutoff.Area(1)) & (x(:,2)>cutoff.Solidity(1));
   fcn.soft_pass = @(x) (x(:,1) < cutoff.Area(2)) & (x(:,1) > cutoff.Area(1)) & (x(:,2)>cutoff.Solidity(2));
end
   fcn.end_pass = @(x) (x(:,1) < cutoff.Area(2)) & (x(:,1) > cutoff.Area(1));

%%
