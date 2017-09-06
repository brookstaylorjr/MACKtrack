function [CellMeasurements, ModuleData] = ppargModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = ppargModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PPARGMODULE measures PPARg intensity (similar to nuclear intensity module, but with modified
% normalization routine. If a 2nd aux image is passed, it will be saved as "CEBP"
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

iteration  = ModuleData.iter;
measure_cc = label2cc(labels.Nucleus,0);

% Normalization 1: flatfield correction -> estimate range for scaling
corr_img = flatfieldcorrect(double(AuxImages{1})-double(parameters.Flatfield{end}),double(parameters.Flatfield{1}));

% Normalization 2: mode-balance - bimodal distribution assumed after dropping nuclei (leaves cytoplasmic + b.g.)
corr_img = corr_img - min(corr_img(:));
tmp = corr_img;
tmp(imdilate(labels.Nucleus>0,diskstrel(round(parameters.MinNucleusRadius*2)))) = []; % Drop foreground objects for correction calculation
[~, dist1] = modebalance(tmp,2,ModuleData.BitDepth,'measure');
corr_img = (corr_img - dist1(1)); % Background subtract (DON'T divide)
AuxImages{1} = corr_img;


% On first call, initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'MeanPPARg')
    % Intensity-based measurement initialization
    CellMeasurements.MeanPPARg =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.IntegratedPPARg =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MedianPPARg = nan(parameters.TotalCells,parameters.TotalImages);
end

% Cycle through each cell and assign measurements
for n = 1:measure_cc.NumObjects
    CellMeasurements.MeanPPARg(n,iteration) = nanmean(AuxImages{1}(measure_cc.PixelIdxList{n}));
    CellMeasurements.IntegratedPPARg(n,iteration) = nansum(AuxImages{1}(measure_cc.PixelIdxList{n}));
    CellMeasurements.MedianPPARg(n,iteration) = nanmedian(AuxImages{1}(measure_cc.PixelIdxList{n}));
end


% Measure nuclei in 2nd auxiliary, if it is specified
if ~isempty(AuxImages{2})
    corr_img = flatfieldcorrect(double(AuxImages{2})-double(parameters.Flatfield{end}),double(parameters.Flatfield{2}));

    % Normalization 2: mode-balance - bimodal distribution assumed after dropping nuclei (leaves cytoplasmic + b.g.)
    corr_img = corr_img - min(corr_img(:));
    tmp = corr_img;
    tmp(imdilate(labels.Nucleus>0,diskstrel(round(parameters.MinNucleusRadius*2)))) = []; % Drop foreground objects for correction calculation
    [~, dist1] = modebalance(tmp,2,ModuleData.BitDepth,'measure'); 
    corr_img = (corr_img - dist1(1)); % Background subtract (DON'T divide)
    AuxImages{2} = corr_img;


    % On first call, initialize all new CellMeasurements fields 
    if ~isfield(CellMeasurements,'MeanCEBP')
        % Intensity-based measurement initialization
        CellMeasurements.MeanCEBP =  nan(parameters.TotalCells,parameters.TotalImages);
        CellMeasurements.IntegratedCEBP =  nan(parameters.TotalCells,parameters.TotalImages);
        CellMeasurements.MedianCEPB = nan(parameters.TotalCells,parameters.TotalImages);
    end

    % Cycle through each cell and assign measurements
    for n = 1:measure_cc.NumObjects
        CellMeasurements.MeanCEBP(n,iteration) = nanmean(AuxImages{2}(measure_cc.PixelIdxList{n}));
        CellMeasurements.IntegratedCEBP(n,iteration) = nansum(AuxImages{2}(measure_cc.PixelIdxList{n}));
        CellMeasurements.MedianCEPB(n,iteration) = nanmedian(AuxImages{2}(measure_cc.PixelIdxList{n}));
    end

end