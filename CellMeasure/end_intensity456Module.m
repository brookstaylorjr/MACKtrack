function [CellMeasurements, ModuleData] = end_intensity456Module(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = end_intensity456Module(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% END_INTENSITYMODULE performs endopoint expression measurement (e.g. post staining of the same cells being imaged)
% Module will try to identify cell boundaries (if not already identified in tracking) from AuxImages{1}, then calls
% intensityModule on this last frame.
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
    
    % If cells were not segmented, use the new cell label matrix identified in end_intensityModule
    if strcmpi(parameters.ImageType,'none')
        labels.Cell = ModuleData.CellLabel_new;
    end

    % Make measurements (make sure to trigger whole-cell measurments as appropriate)
    tmp_measurements = struct; 
    tmp_params = parameters; tmp_params.TotalImages = 1; tmp_params.ImageType = 'fluorescence';
    tmp_data = ModuleData; tmp_data.iter = 1;
    tmp_measurements = intensity456Module(tmp_measurements,tmp_params, labels, AuxImages, tmp_data);  
    names = fieldnames(tmp_measurements); 
    
    for m = 1:length(names)
        CellMeasurements.(['end_',names{m}]) = tmp_measurements.(names{m})(:,end);
    end
 
end
    
