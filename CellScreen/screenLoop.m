function [] = screenLoop(parameters)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [] = screenLoop(parameters)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SCREENLOOP crawls through all subdirectories specified under parameters.ImagePath, looking for layout.xlsx files that
% correspond to plate screener images. Cells in correponding channels are segmented and measured. The input file 
% structure will be echoed under the output (save) directory.
%
% (See the subfunction 'wellmatch' for valid naming schema)
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Get locations.mat, identify starting directory
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
load([home_folder(1:slash_idx(end-1)), 'locations.mat'],'-mat')

start_dir = namecheck([locations.scope, filesep,parameters.ImagePath]);

% Crawl down from starting directory - find all plate layouts (and corresponding image folders)
[layout_dir, image_dir] = getImageDirectories(start_dir);
save_dir = [locations.data,filesep,parameters.SaveDirectory];
mkdir(save_dir)

% Pre-save parameters for each tracked set
for idx = 1:length(layout_dir)
    save_subdir = [save_dir,filesep,layout_dir{idx}(length(start_dir)+1:end)];
    mkdir(save_subdir)
    parameters.ImagePath_full = namecheck(image_dir{idx}(length(locations.scope)+1:end));
    save([save_subdir,filesep,'TrackingParameters.mat'],'parameters')
end

% Process flatfields - stash original values
params_save = parameters;
if isfield(parameters,'Flatfield')
    params_save.Flatfield_orig = parameters.Flatfield;
    parameters.Flatfield = processFlatfields(parameters.Flatfield);
end


for idx = 1:length(layout_dir)
    % a) PER PLATE LAYOUT: identify al conditions (and corresponding wells) to be measured
    save_subdir = [save_dir,filesep,layout_dir{idx}(length(start_dir)+1:end)];
    mkdir([save_subdir,filesep,'SegmentedImages'])
    mkdir([save_subdir,filesep,'NuclearLabels'])
    if ~strcmpi(parameters.ImageType,'none')
            mkdir([save_subdir,filesep,'CellLabels'])
    end
    AllData = struct;
    [conditions, wells] = parseLayout(layout_dir{idx});
    
    % Get image list, and associated bit depth of full resolution (i.e. not thumbnail) images
    image_names = quickdir(image_dir{idx});
    [tmp_names, parameters.scope_type] = wellmatch(image_names, wells{1}{1});
    if isempty(tmp_names)
        error(['No images found that correspond to well ',wells{1}{1},'. Double-check "layout.xlsx" file'])
    end
    
    % Get image bit depth (for processing)
    imfo = imfinfo([image_dir{idx}, filesep, tmp_names{1}]);
    parameters.BitDepth = imfo.BitDepth;
    
    
    % Analyze all wells/images per conditon (in parallel, if selected)
    Data = cell(size(conditions));
    if parameters.Parallel
        parfor k = 1:length(conditions)
            Data{k} = screenprocess(image_dir{idx}, image_names, wells{k}, save_subdir, parameters);
        end
    else
        for k = 1:length(conditions)
            Data{k} = screenprocess(image_dir{idx}, image_names, wells{k}, save_subdir, parameters);
        end
    end
    % Combine data into a more "useful" structure
    for k = 1:length(conditions)
        condition_field = matlab.lang.makeValidName(conditions{k});
        AllData.(condition_field) = Data{k};
    end
    
    % Save full measurments file
    save([save_subdir,filesep,'AllData.mat'],'AllData')

    % Reformat data in a flattened format for R/Python users as well
    [AllData_mat, AllData_fields, AllData_conditions] = AllData_to_mat(AllData);
    save([save_subdir,filesep,'AllData_R.mat'],'AllData_mat','AllData_fields','AllData_conditions')

end

        