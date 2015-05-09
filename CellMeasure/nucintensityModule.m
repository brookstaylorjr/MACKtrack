function [CellMeasurements, ModuleData] = nucintensityModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% INTENSITYMODULE measures cellular intensity in AuxImages{1} channel
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

iteration  = ModuleData.iter;

% Mode-balance 1st auxililiary image - unimodal distribution assumed
    if ~isfield(ModuleData,'distr')
        [~, ModuleData.distr] = modebalance(AuxImages{1},1,ModuleData.BitDepth,'measure'); 
    else
        AuxImages{1} = modebalance(AuxImages{1},1,ModuleData.BitDepth,'correct',ModuleData.distr);
    end


% On first call, initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'MeanIntensityNuc')
    % Intensity-based measurement initialization
    CellMeasurements.MeanIntensityNuc =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.IntegratedIntensityNuc =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MedianIntensityNuc = nan(parameters.TotalCells,parameters.TotalImages);
end

% Cycle through each cell and assign measurements
cells = unique(labels.Nucleus(labels.Nucleus>0));
for n = 1:length(cells)
    CellMeasurements.MeanIntensityNuc(cells(n),iteration) = nanmean(AuxImages{1}(labels.Nucleus==cells(n)));
    CellMeasurements.IntegratedIntensityNuc(cells(n),iteration) = nansum(AuxImages{1}(labels.Nucleus==cells(n)));
    CellMeasurements.MedianIntensityNuc(cells(n),iteration) = nanmedian(AuxImages{1}(labels.Nucleus==cells(n)));
end

ModuleDataOut = ModuleData;