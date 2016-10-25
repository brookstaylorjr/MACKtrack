function [CellMeasurements, ModuleData] = ppargModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% [CellMeasurements, ModuleData] = ppargModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PPARGMODULE measures PPARg intensity (similar to nuclear intensity module, but with modified
% normalization routine.
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

iteration  = ModuleData.iter;
% Mode-balance 1st auxililiary image - trimodal distribution assumed (nuclear expression, cytoplasmic expression, and b.g.)
if ~isfield(ModuleData,'distr')
    [~, ModuleData.distr] = modebalance(AuxImages{1},3,ModuleData.BitDepth,'measure'); 
end
% Background subtract - don't divide to correct for standard deviation change
% (Confluency makes this measurement slightly unreliable).
AuxImages{1} = AuxImages{1} - ModuleData.distr(1);

% On first call, initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'MeanIntensityNuc')
    % Intensity-based measurement initialization
    CellMeasurements.MeanPPARg =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.IntegratedPPARg =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MedianPPARg = nan(parameters.TotalCells,parameters.TotalImages);
end

% Cycle through each cell and assign measurements
cells = unique(labels.Nucleus(labels.Nucleus>0));
for n = 1:length(cells)
    CellMeasurements.MeanPPARg(cells(n),iteration) = nanmean(AuxImages{1}(labels.Nucleus==cells(n)));
    CellMeasurements.IntegratedPPARg(cells(n),iteration) = nansum(AuxImages{1}(labels.Nucleus==cells(n)));
    CellMeasurements.MedianPPARg(cells(n),iteration) = nanmedian(AuxImages{1}(labels.Nucleus==cells(n)));
end

ModuleDataOut = ModuleData;