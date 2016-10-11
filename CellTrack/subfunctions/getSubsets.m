function [subset_list] = getSubsets(obj, all_neighbors, n)

neighb1 = filterNeighbors(all_neighbors{obj},obj);

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
        merged_neighbors = filterNeighbors(merged_neighbors,check_subset{1});
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


function filtered_set = filterNeighbors(set,self)
set(ismember(set,self)) = [];
set(set==0) = [];
filtered_set = unique(set);
