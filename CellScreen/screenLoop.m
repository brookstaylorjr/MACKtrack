function [] = screenLoop(parameters)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [] = screenLoop(parameters)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SCREENLOOP crawls through all subdirectories specified under ImagePath, looking for layout.xlsx files that correspond
% to plate screener images. Cells in correponding channels are segmented and measured. The input file structure will
% be echoed under the output (save) directory.
%
% Function is currently configured to process default naming from a MicroXLS . 
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Get locations.mat, identify starting directory
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
load([home_folder(1:slash_idx(end-1)), 'locations.mat'],'-mat')

start_dir = [locations.scope, filesep,parameters.ImagePath];

% Crawl down from starting directory - find all plate layouts (and corresponding image folders)
[layout_dir, image_dir, all_dir] = getImageDirectories(start_dir);
save_dir = [locations.data,filesep,parameters.SaveDirectory];
mkdir(save_dir)


for idx = 1:length(layout_dir)
    % a) PER PLATE LAYOUT: identify al conditions (and corresponding wells) to be measured
    disp(['Analyzing subdirectory "',layout_dir{1}(length(start_dir)+1:end),'"...'])
    save_subdir = [save_dir,filesep,layout_dir{idx}(length(start_dir)+1:end)];
    mkdir(save_subdir)
    mkdir([save_subdir,filesep,'SegmentedImages'])
    mkdir([save_subdir,filesep,'NuclearLabels'])
    AllData = struct;
    [conditions, wells] = parseLayout(layout_dir{idx});
    % Get image list, and associated bit depth of full resolution (i.e. not thumbnail) images
    image_names = struct2cell(dir(image_dir{idx}));
    image_names = image_names(1,:)';
    tmp_name = image_names{find(cellfun(@isempty,strfind(image_names,'thumb'))...
        &~cellfun(@isempty,strfind(image_names,['_',wells{1}{1},'_'])),1,'first')};
    imfo = imfinfo([image_dir{idx}, filesep, tmp_name]);
    parameters.BitDepth = imfo.BitDepth;
    
    % Analyze all wells/images per conditon (in parallel)
    Data = cell(size(conditions));
    parfor k = 1:length(conditions)
        Data{k} = microxlprocess(image_dir{idx}, image_names, wells{k}, save_subdir, parameters);
    end
    % Combine data into a more "useful" structure
    for k = 1:length(conditions)
        condition_field = matlab.lang.makeValidName(conditions{k});
        AllData.(condition_field) = Data{k};
    end
    
    % Save full measurments file
    save([save_subdir,filesep,'AllData.mat'],'AllData')
end

        