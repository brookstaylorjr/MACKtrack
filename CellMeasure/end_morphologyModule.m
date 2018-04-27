function [CellMeasurements, ModuleData] = end_morphologyModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% end_morphologyModule measures simple aspects of nuclei/cellular
% morphology using only the last frame
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
    
    % If cells were not segmented, use 1st available aux image to identify cell boundaries
    if strcmpi(parameters.ImageType,'none')
        idx = 1;
        aux_image = AuxImages{idx};
        while isempty(aux_image)
            idx = idx+1;
            aux_image = AuxImages{idx};
            if idx>10
                break;
            end
        end
        parameters.CellFF = idx; % Default to (corresponding) flatfield
        data = fluorescenceID(aux_image, parameters, []);
        data.nuclei = labels.Nucleus;
        tmp_out = fluorescenceSegment(data, aux_image, parameters);
        labels.Cell = tmp_out.cells;
    end
    
     % Make measurements (make sure to trigger whole-cell measurments as appropriate)
    tmp_measurements = struct; 
    tmp_params = parameters; tmp_params.TotalImages = 1; tmp_params.ImageType = 'fluorescence';
    tmp_data = ModuleData; tmp_data.iter = 1;
    AuxImages{1} = AuxImages{find(~isempty(AuxImages),1,'last')};
    
    tmp_measurements = morphologyModule(tmp_measurements,tmp_params, labels, AuxImages, tmp_data);  
    names = fieldnames(tmp_measurements); 
    
    for m = 1:length(names)
        CellMeasurements.(['end_',names{m}]) = tmp_measurements.(names{m})(:,end);
    end
end