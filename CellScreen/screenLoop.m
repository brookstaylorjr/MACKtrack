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

start_dir = ffp([locations.scope, filesep,parameters.ImagePath]);

% Crawl down from starting directory - find all plate layouts (and corresponding image folders)
[layout_dir, image_dir, all_dir] = getImageDirectories(start_dir);
save_dir = [locations.data,filesep,parameters.SaveDirectory];
mkdir(save_dir)


for idx = 1:length(layout_dir)
    % a) PER PLATE LAYOUT: identify al conditions (and corresponding wells) to be measured
    save_subdir = [save_dir,filesep,layout_dir{idx}(length(start_dir)+1:end)];
    mkdir(save_subdir)
    mkdir([save_subdir,filesep,'SegmentedImages'])
    mkdir([save_subdir,filesep,'NuclearLabels'])
    if ~strcmpi(parameters.ImageType,'none')
            mkdir([save_subdir,filesep,'CellLabels'])
    end
    AllData = struct;
    [conditions, wells] = parseLayout(layout_dir{idx});
    % Get image list, and associated bit depth of full resolution (i.e. not thumbnail) images
    image_names = quickdir(image_dir{idx});
    try
    tmp_name = image_names{find(cellfun(@isempty,strfind(image_names,'thumb'))...
        &~cellfun(@isempty,strfind(image_names,['_',wells{1}{1},'_'])),1,'first')};
    catch me
        error(['No images found that correspond to well ',wells{1}{1},'. Double-check "layout.xlsx" file'])
    end
        
    % Make parameter additions: 1) bit depth and 2) flatfield image fits.
    imfo = imfinfo([image_dir{idx}, filesep, tmp_name]);
    parameters.BitDepth = imfo.BitDepth;
    % Calculate flatfield images - replace them in parameters
    if isfield(parameters,'Flatfield')
        X = backgroundcalculate([imfo.Height,imfo.Width]);
        warning off MATLAB:nearlySingularMatrix
        for i = 1:length(parameters.Flatfield)
            corr_img = parameters.Flatfield{i};
            pStar = (X'*X)\(X')*corr_img(:);
            % Apply correction
            corr_img = reshape(X*pStar,size(corr_img));
            parameters.Flatfield{i} = corr_img-min(corr_img(:));
        end
    end
    
    % Analyze all wells/images per conditon (in parallel, if selected)
    Data = cell(size(conditions));
    if parameters.Parallel
        parfor k = 1:length(conditions)
            Data{k} = microxlprocess(image_dir{idx}, image_names, wells{k}, save_subdir, parameters);
        end
    else
        for k = 1:length(conditions)
            Data{k} = microxlprocess(image_dir{idx}, image_names, wells{k}, save_subdir, parameters);
        end
    end
    % Combine data into a more "useful" structure
    for k = 1:length(conditions)
        condition_field = matlab.lang.makeValidName(conditions{k});
        AllData.(condition_field) = Data{k};
    end
    
    % Save full measurments file
    save([save_subdir,filesep,'AllData.mat'],'AllData')
    % Save tracking parameters file - adjust to include full image directory
    parameters.ImagePath_full = image_dir{idx};
    save([save_subdir,filesep,'TrackingParameters.mat'],'parameters')
    
    % Reformat data in a flattened format for R/Python users as well
    [AllData_mat, AllData_fields, AllData_conditions] = AllData_to_mat(AllData);
    save([save_subdir,filesep,'AllData_R.mat'],'AllData_mat','AllData_fields','AllData_conditions')

end

        