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

% Mode-balance 1st auxililiary image - trimodal distribution assumed (nuclear expression, cytoplasmic expression, and b.g.)
[~, dist1] = modebalance(AuxImages{1},3,ModuleData.BitDepth,'measure'); 
% Background subtract - don't divide to correct for standard deviation change
% (Confluency makes this stdev measurement slightly unreliable).
AuxImages{1} = AuxImages{1} - dist1(1);

% On first call, initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'MeanIntensityNuc')
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
    % Mode-balance 1st auxililiary image - trimodal distribution assumed (nuclear expression, cytoplasmic expression, and b.g.)
    [~, dist2] = modebalance(AuxImages{2},3,ModuleData.BitDepth,'measure'); 
    % Background subtract - don't divide to correct for standard deviation change
    % (Confluency makes this stdev measurement slightly unreliable).
    AuxImages{2} = AuxImages{2} - dist2(1);

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