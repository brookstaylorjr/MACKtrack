function [CellMeasurements, ModuleData] = penult_intensityModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = penult_intensityModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PENULT_INTENSITYMODULE performs expression measurement (e.g. on post staining of the same cells being 
% imaged). These images are assumed to be in the PENULTIMATE timepoint, but cell mask will be found in FINAL
% timepoint. (This allows tracking through a set of post-stained/re-stained images)
%
%
% Though it refers to a timepoint previous, it triggers AFTER end_intensityModule, since it relies on a 
% label matrix made by 
% 
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


%% [Engage on last frame only]
if ModuleData.iter == parameters.TotalImages
    
    % This module assumes cells were not segmented in tracking, use the cell label matrix from end_intensityModule
    labels.Cell = ModuleData.CellLabel_new;
    
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
    
    % Make measurements (make sure to trigger whole-cell measurements as appropriate)
    tmp_measurements = struct; 
    tmp_params = parameters; tmp_params.TotalImages = 1; tmp_params.ImageType = 'fluorescence';
    tmp_data = ModuleData; tmp_data.iter = 1;
    tmp_measurements = intensityModule(tmp_measurements,tmp_params, labels_new, AuxImages, tmp_data);  
    names = fieldnames(tmp_measurements); 
    
    for m = 1:length(names)
        CellMeasurements.(['penult_',names{m}]) = tmp_measurements.(names{m})(:,end);
    end
 
end
    
