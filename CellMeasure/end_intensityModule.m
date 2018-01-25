function [CellMeasurements, ModuleData] = end_intensityModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = end_intensityModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
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
        
        % Save a diagnostic output version of this image
        home_folder = mfilename('fullpath');
        slash_idx = strfind(home_folder,filesep);
        load([home_folder(1:slash_idx(end-1)), 'locations.mat'],'-mat')
        save_dir = namecheck([locations.data,filesep,parameters.SaveDirectory,filesep,'EndpointSegmentation',filesep]);
        if ~exist(save_dir,'dir');  mkdir(save_dir); end
        tmp1 = aux_image;
        tmp1(tmp1==min(tmp1(:))) = [];
        tmp1(tmp1==max(tmp1(:))) = [];
        tmp1 = modebalance(tmp1,0, ModuleData.BitDepth,'display');               
        % Non-confluent case - set low saturation @ 3xS.D. below bg level
        if parameters.Confluence ~= 1

            pct = 90:.5:99;
            hi_val = prctile(tmp1,pct);       
            saturation_val = [-3 prctile(tmp1,1+findelbow(pct,hi_val))];
            alpha = 0.4;
        else % Confluent case: unimodal distribution is foreground - use a different lower limit.
            saturation_val = [-4 prctile(tmp1(:),90)];
            alpha = 0.55;
        end
        saveFig(aux_image,labels.Cell,labels.Nucleus,[],ModuleData.BitDepth,...
            [save_dir,'Endpoint-pos_',numseq(ModuleData.i,3),'.jpg'],alpha,[1024 1024], saturation_val);
    end

    % Make measurements (make sure to trigger whole-cell measurments as appropriate)
    tmp_measurements = struct; 
    tmp_params = parameters; tmp_params.TotalImages = 1; tmp_params.ImageType = 'fluorescence';
    tmp_data = ModuleData; tmp_data.iter = 1;
    tmp_measurements = intensityModule(tmp_measurements,tmp_params, labels, AuxImages, tmp_data);  
    names = fieldnames(tmp_measurements); 
    
    for m = 1:length(names)
        CellMeasurements.(['end_',names{m}]) = tmp_measurements.(names{m})(:,end);
    end
 
end
    
