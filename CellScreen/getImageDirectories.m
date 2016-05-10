function [layout_dir, image_dir, all_dir] = getImageDirectories(start_dir)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%[layout_dir, image_dir, all_dir] = getImageDirectories(start_dir)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% GETIMAGEDIRECTORIES sweeps through all subdirectories of a specified start_dir and identifies all screening
% experiments (should contain a 'layout.xlsx' file), as well as the associated image directory with that plate's 
% layout.
%
% INPUTS
% start_dir  starting directory - a single structure will be made, containing fields & subfields corresponding
%            to the file hierarchy underneath
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%% Iterate through folder structure: look for layout. If not found, look for subfolders and work downwards 
layout_check  = @(dir1) exist([dir1,filesep,'layout.xlsx'],'file');
filter_dir = @(dir_cell) find((cell2mat(dir_cell(4,:))) & ~strcmp(dir_cell(1,:),'.') & ~strcmp(dir_cell(1,:),'..'));
found_layout = layout_check(start_dir);
all_dir = {start_dir};
layout_dir = {};

% Search ALL subdirectories and add path (until there are no more directories to search)
search_dir = all_dir;
while length(search_dir)>0
    curr_dir = search_dir{1};
    % Look for layout.xlsx - if it's found, add to layout_dir (and don't continue looking for subdirectories there)
    if layout_check(curr_dir)
        layout_dir = cat(1,layout_dir,curr_dir);
    else
        % layout.xlsx not found: find/add all subdirectories to be searched
        dir_cell = struct2cell(dir(curr_dir));
        dir_idx = filter_dir(dir_cell);
        % Add all directories to be searched
        for i = 1:length(dir_idx)
            all_dir = cat(1,all_dir,[curr_dir,filesep,dir_cell{1,dir_idx(i)}]);
            search_dir = cat(1,search_dir,[curr_dir,filesep,dir_cell{1,dir_idx(i)}]);
        end
    end
    % Knock out current search directory
    search_dir(1) = [];
end

% For each layout directory, find the corresponding image directory (i.e. first directory we find that has .tif images)
image_dir = cell(size(layout_dir));
for i = 1:length(layout_dir)
    search_dir = {layout_dir{i}};
    while length(search_dir)>0
        curr_dir = search_dir{1};
        dir_cell = struct2cell(dir(curr_dir));
        % Look for (at least two) .tif files. If it's found, add it to image_dir (and break)
        if sum(~cellfun(@isempty,strfind(dir_cell(1,:),'.tif')))>1
            image_dir{i} = curr_dir;
            break
        else % No .tif files - keep searching
            dir_idx = filter_dir(dir_cell);
            % Add all directories to be searched
            for j = 1:length(dir_idx)
                search_dir = cat(1,search_dir,[curr_dir,filesep,dir_cell{1,dir_idx(j)}]);
            end
        end
        search_dir(1) = [];
    end

end

