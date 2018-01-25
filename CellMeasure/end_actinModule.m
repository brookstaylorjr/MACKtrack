function [CellMeasurements, ModuleData] = end_actinModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = end_actinModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% END_ACTINMODULE identifies actin fibers in labeled cells, quantifying their orientation and scoring a cell by the 
% relativel level of organized fibers.
%
% Module will try to identify cell boundaries (if not already identified in tracking) from AuxImages{1}, then calls
% actinModule on this last frame. If 2 images are provided, it will be assumed that the LAST one is the actin image. 
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
        ModuleData.CellLabel_new = labels.Cell;
    end


    % Set up a save directory for diagnostic images
    home_folder = mfilename('fullpath');
    slash_idx = strfind(home_folder,filesep);
    load([home_folder(1:slash_idx(end-1)), 'locations.mat'],'-mat')
    save_dir = namecheck([locations.data,filesep,parameters.SaveDirectory,filesep,'End_ActinID',filesep]);
    if ~exist(save_dir,'dir');  mkdir(save_dir); end

    
    
    % Make measurements (make sure to trigger whole-cell measurments as appropriate)
    tmp_measurements = struct; 
    tmp_params = parameters; tmp_params.TotalImages = 1; tmp_params.ImageType = 'fluorescence';
    tmp_data = ModuleData; tmp_data.iter = 1; tmp_data.save_dir = save_dir;
    AuxImages{1} = AuxImages{find(~isempty(AuxImages),1,'last')};
    
    tmp_measurements = actinModule(tmp_measurements,tmp_params, labels, AuxImages, tmp_data);  
    names = fieldnames(tmp_measurements); 
    
    for m = 1:length(names)
        CellMeasurements.(['end_',names{m}]) = tmp_measurements.(names{m})(:,end);
    end
    
    
    
 
    
    
end
    
