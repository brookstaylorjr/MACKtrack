function [CellMeasurements, ModuleData] = multi_fretModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = multi_fretModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MULTI_FRETMODULE measures a FRET ratio (after correction) between FRET and CFP images. Additionally, it will scan for
% subseampled images (i.e. multiple measurement image pairs per tracking set)
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


% On first call, initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'MeanFRET_multi')
    % Intensity-based measurement initialization
    CellMeasurements.MeanFRET_multi =  nan(parameters.TotalCells,5*parameters.TotalImages);
    CellMeasurements.IntegratedFRET_multi =  nan(parameters.TotalCells,5*parameters.TotalImages);
    CellMeasurements.MedianFRET_multi = nan(parameters.TotalCells,5*parameters.TotalImages);
    ModuleData.col_idx = 0;
end

% Skip last image
if ModuleData.iter == parameters.TotalImages
    disp('Skipping last image.')
    CellMeasurements.MeanFRET_multi(:,ModuleData.col_idx+1:end) = [];
    CellMeasurements.IntegratedFRET_multi(:,ModuleData.col_idx+1:end) = [];
    CellMeasurements.MedianFRET_multi(:,ModuleData.col_idx+1:end) = [];
    return;
end


% Define label matricies for cells/nuclei
tmp_label = labels.Cell;
tmp_label(labels.Nucleus>0) = 0; % cytoplasm only
measure_cc = label2cc(tmp_label,0);




% Get total number of image pairs we're going to measure - scan for images with same timepoints
% i.e. if we're at t100, look for existence of t100.01, t100.02, etc.

t_expr = '[0-9]_T[0-9]*_'; % regexp will match this in finding t vals -> note: value starts 3 pos past the output idx

name1 = ModuleData.AuxName{1};
name2 = ModuleData.AuxName{2};
[~,idx1] = regexp(name1,t_expr);
[~,idx2] = regexp(name2,t_expr);

for j = 1:5
    test_name1 = [name1(1:idx1-1),'.',numseq(j,2),name1(idx1:end)];
    test_name2 = [name2(1:idx2-1),'.',numseq(j,2),name2(idx2:end)];
    if exist(test_name1,'file') && exist(test_name2,'file')
        AuxImages{1} = cat(3,AuxImages{1},checkread(test_name1,ModuleData.BitDepth));
        AuxImages{2} = cat(3,AuxImages{2},checkread(test_name2,ModuleData.BitDepth));
    end
end


% Get ratiometric measurement for FRET
for i = 1:size(AuxImages{1},3)
    ModuleData.col_idx = ModuleData.col_idx+1;
    fret = AuxImages{1}(:,:,i); 
    fret = fret - double(parameters.Flatfield{end});
    fret = flatfieldcorrect(fret,double(parameters.Flatfield{1}));
    fret = fret-prctile(fret(:),2); % Background subtract
    fret(fret<16) = 1; % add floor to image 
    
    cfp = AuxImages{2}(:,:,i);
    cfp = cfp - double(parameters.Flatfield{end});
    cfp = flatfieldcorrect(cfp,double(parameters.Flatfield{1}));
    cfp(cfp<16) = 16; % add floor to image 
    fret_image = (fret)./(cfp);

    % Cycle through each cell and assign measurements
    for n = 1:measure_cc.NumObjects
        CellMeasurements.MeanFRET_multi(n,ModuleData.col_idx) = nanmean(fret_image(measure_cc.PixelIdxList{n}));
        CellMeasurements.IntegratedFRET_multi(n,ModuleData.col_idx) = nansum(fret_image(measure_cc.PixelIdxList{n}));
        CellMeasurements.MedianFRET_multi(n,ModuleData.col_idx) = nanmedian(fret_image(measure_cc.PixelIdxList{n}));
    end
end
