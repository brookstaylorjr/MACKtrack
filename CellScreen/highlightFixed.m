function [] = highlightFixed(CellData, image_list, track_folder)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [] = highlightFixed(CellData, image_list, track_folder)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% HIGHLIGHTFIXED will highlight a (provided) subset of cells, showing their nuclear (and cell, if applicable)  
% boundaries on the originally-aquired images. If provided subset spans several sites, a new figure window will
% be opened for each site.
%
% INPUTS
% CellData      selected cells, from AllData.(condition).CellData. Contains: [ well idx | img_idx | cell_idx ]
% image_list    segmented image list, from AllData.(condition).Images
% track_folder  path to corresponding analysis (i.e. the folder that contains 'AllData.mat' and 
%                    'TrackingParameters.mat'
%
%%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



%%

% Load scope/data locations
home_folder = mfilename('fullpath');
slash_idx = strfind(home_folder,filesep);
home_folder = home_folder(1:slash_idx(end-1));
load([home_folder, 'locations.mat'],'-mat')
if ~exist(locations.scope)
    error(['Invalid mount location for images: "',locations.scope,...
        '" not found. Please load, update, and re-save "locations.mat"'])
end


% Load parameters
disp(['Loading in "TrackingParameters.mat" from ', track_folder])
load(namecheck([track_folder,filesep,'TrackingParameters.mat']))


% Get full contents of image directory
handles.image_dir = namecheck([locations.scope, filesep, parameters.ImagePath_full]);
tif_images = quickdir(handles.image_dir);
tif_images(cellfun(@isempty,strfind(tif_images,'.tif'))) = []; % Drop non-tif entries

% Scrape well/position from supplied image list. Open up a new fig for each frame that's included
image_names = image_list(unique(CellData(:,2)),1);
image_idx = unique(CellData(:,2));



% Get all fluorescent channels 


for i = 1:length(image_names)
    cell_idx = CellData(CellData(:,2)==image_idx(i),3);
    idx = strfind(image_names{i},'_');
    matchstr = image_names{i}(1:idx(end));
    match_names = tif_images(~cellfun(@isempty,strfind(tif_images,matchstr)));
    match_names = match_names(cellfun(@isempty,strfind(match_names,'_thumb'))); % Drop thumbnails
    
    % Load corresponding nuclear label
    well = image_names{i}(1+idx(end-2):idx(end)-1);
    well = [well(1:3),'_',numseq(eval(well(6:end)),2)];
    load(namecheck([track_folder,filesep,'NuclearLabels',filesep,'NuclearLabel-',well,'.mat'])) 
    NuclearLabel(~ismember(NuclearLabel, cell_idx)) = 0; % Drop all except selected cells
    
    % Load /filter CellLabel, if it exists
    add_cell = 0;
    name1 = namecheck([track_folder,filesep,'CellLabels',filesep,'CellLabel-',well,'.mat']);
    if exist(name1,'file')
        load(name1);
        CellLabel(~ismember(CellLabel, cell_idx)) = 0; % Drop all except selected cells
        add_cell = 1;
    end
    
    
    % Make Z stack with all channels.
    
    
    z = [];
    
    for j = 1:length(match_names)
        img = checkread(namecheck([locations.scope,filesep,parameters.ImagePath_full,filesep,match_names{j}]));
        cap = prctile(img(:),[2 98]);
        img = (img-min(cap))/range(cap);
        img(img<0) = 0;
        img(img>1) = 1;  
        img((imdilate(NuclearLabel,ones(3))-NuclearLabel)>0) = max(img(:));
        if add_cell
            img((imdilate(CellLabel,ones(3))-CellLabel)>0) = 0.7*max(img(:));
        end
        z = cat(3,z,img);   
    end
    colormaps = loadcolormaps;
    figure,multistack(z), colormap(colormaps.purple_hot)
    
    
    
end






