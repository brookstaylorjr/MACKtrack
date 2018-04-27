function [CellMeasurements, ModuleData] = penult_actinModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = penult_actinModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PENULT_ACTINMODULE identifies actin fibers in labeled cells, quantifying their orientation and scoring a cell by the 
% relativel level of organized fibers. These images are assumed to be in the PENULTIMATE timepoint, but cell mask will 
% be found in FINAL timepoint. (This allows tracking through a set of post-stained/re-stained images, that may be 
% shifted slightly)
%
% Though the measurement is made on the 2nd-to-last set of images, it isn't triggered until AFTER end_intensityModule, 
% since it relies on a label matrix made by this module -> end_intensityModule MUST be selected.
%
% OUTPUTS:
% CellMeasurements    structure with fields corresponding to cell measurements
%
% INPUTS:
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


%% [Engage on last frame only]
if ModuleData.iter == parameters.TotalImages
    
    % If cells were not segmented, use the new cell label matrix identified in end_intensityModule
    if strcmpi(parameters.ImageType,'none')
        labels.Cell = ModuleData.CellLabel_new;
    end
    
    % (Likewise, the Nuclear label matrix will be the one corresponding to the final/ground-truth measurement
    
    % Load offsets. Apply to masks existing cell/nuclear label matricies
    load([parameters.XYDir,filesep,'ImageJumps.mat'])
    offset_post = image_jumps(image_jumps(:,1)== parameters.TimeRange(ModuleData.iter),2:3);
    offset_pre = image_jumps(image_jumps(:,1)== parameters.TimeRange(ModuleData.iter-1),2:3);
    if isempty(offset_post) || isempty(offset_pre)
        error('Image jumps weren''t calculated for post-stain points - retracking is required')
    end 
    adj = (offset_pre-offset_post);
    
    adj_labels = {labels.Nucleus,labels.Cell};
    for idx = 1:length(adj_labels)
        tmp_label = adj_labels{idx};
        % Adjust rows
        if adj(1) > 0
            tmp_label = [zeros(adj(1),size(tmp_label,2)); tmp_label(1:end-adj(1),:) ];    
        elseif adj(1) < 0
            tmp_label = [tmp_label(-adj(1)+1:end,:), zeros(adj(1),size(tmp_label,2))];
        end
        % Adjust cols
        if adj(2) > 0
            tmp_label = [zeros(size(tmp_label,1),adj(2)), tmp_label(:,1:end-adj(2)) ];    
        elseif adj(2) < 0
            tmp_label = [tmp_label(:,-adj(2)+1:end), zeros(size(tmp_label,1),-adj(2))];
        end
        adj_labels{idx} = tmp_label;
    end
    
    labels_new.Nucleus = adj_labels{1};
    labels_new.Cell = adj_labels{2};
    
    % Set up a save directory for diagnostic images
    home_folder = mfilename('fullpath');
    slash_idx = strfind(home_folder,filesep);
    load([home_folder(1:slash_idx(end-1)), 'locations.mat'],'-mat')
    save_dir = namecheck([locations.data,filesep,parameters.SaveDirectory,filesep,'Penult_ActinID',filesep]);
    if ~exist(save_dir,'dir');  mkdir(save_dir); end
    
    
    % Make measurements (make sure to trigger whole-cell measurements as appropriate)
    tmp_measurements = struct; 
    tmp_params = parameters; tmp_params.TotalImages = 1; tmp_params.ImageType = 'fluorescence';
    tmp_data = ModuleData; tmp_data.iter = 1; tmp_data.save_dir = save_dir;
    tmp_measurements = actinModule(tmp_measurements,tmp_params, labels_new, AuxImages, tmp_data);  
    names = fieldnames(tmp_measurements); 
    
    for m = 1:length(names)
        CellMeasurements.(['penult_',names{m}]) = tmp_measurements.(names{m})(:,end);
    end
 
end
    
