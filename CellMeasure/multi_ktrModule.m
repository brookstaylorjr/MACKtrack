function [CellMeasurements, ModuleData] = multi_ktrModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = multi_ktrModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MULTI_KTRMODULE measures a KTR ratio (after correction) between FRET and CFP images. Additionally, it will scan for
% subseampled images (i.e. multiple measurement image pairs per tracking set)
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%% STEP 1: FIND All IMAGES (BETWEEN CURRENT AND NEXT "WHOLE" TIMEPOINT)
% Get total number of image pairs we're going to measure - scan for images with same timepoints
% i.e. if we're at t100, look for existence of t100.01, t100.02, etc.
t_expr = '[0-9]_T[0-9]*_'; % regexp will match this in finding t vals -> note: value starts 3 pos past the output idx
aux_tmp = cell(size(AuxImages));
flag  = 1; m =1;
while flag
    flag = 0;
    for n = 1:length(AuxImages)
        name_tmp = ModuleData.AuxName{n};      
        if ~isempty(name_tmp)
            [ ~,idx1] = regexp(name_tmp,t_expr);
            test_name1 = [name_tmp(1:idx1-1),'.',numseq(m,2),name_tmp(idx1:end)];
            if exist(test_name1,'file') 
                aux_tmp{n} = checkread(test_name1,ModuleData.BitDepth);
                flag = 1;
            end
        end
    end
    m = m+1;
    if flag; AuxImages = cat(1,AuxImages,aux_tmp); end
end

% Add corresponding time sequences for "multi" sequences (assume even spacing until next timepoint)
new_ts = linspace(ModuleData.iter, ModuleData.iter+1,size(AuxImages,1)+1);

if ~isfield(CellMeasurements,'MultiKTR_t'); CellMeasurements.MultiKTR_t = []; end
CellMeasurements.MultiKTR_t = cat(2, CellMeasurements.MultiKTR_t, new_ts(1:end-1));


%% STEP 2: MEASURE ALL IMAGES - CONCATENATE MEASUREMENT
tmp_measurements = struct;
tmp_params = parameters; 
tmp_data = ModuleData; tmp_data.iter  = 1; 
tmp_params.TotalImages = size(AuxImages,1); 
for n = 1:size(AuxImages,1)
    tmp_measurements = ktrModule(tmp_measurements,tmp_params, labels, AuxImages(n,:), tmp_data);
    tmp_data.iter = tmp_data.iter+1;
end

names = fieldnames(tmp_measurements);
if ModuleData.iter==1
    CellMeasurements = combinestructures(CellMeasurements,tmp_measurements);
else
    for m = 1:length(names)
        CellMeasurements.(names{m}) = cat(2,CellMeasurements.(names{m}),tmp_measurements.(names{m}));
    end
end

