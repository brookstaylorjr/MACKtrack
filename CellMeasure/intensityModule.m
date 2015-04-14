function [CellMeasurements, ModuleData] = intensityModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
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

% Mode-balance 1st auxililiary image - bimodal distribution assumed
    if ~isfield(ModuleData,'distr')
        [AuxImages{1}, ModuleData.distr] = modebalance(AuxImages{1},2,ModuleData.BitDepth,'measure'); 
    else
        AuxImages{1} = modebalance(AuxImages{1},2,ModuleData.BitDepth,'correct',ModuleData.distr);
    end


% Initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'MeanIntensity')
    % Intensity-based measurement initialization
    CellMeasurements.MeanIntensity =  zeros(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.IntegratedIntensity =  zeros(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MedianIntensity = zeros(parameters.TotalCells,parameters.TotalImages);
end

% Cycle through each cell and assign measurements
cells = unique(labels.Cell(labels.Cell>0));
for n = 1:length(cells)
    CellMeasurements.MeanIntensity(cells(n),iteration) = mean(AuxImages{1}(labels.Cell==cells(n)));
    CellMeasurements.IntegratedIntensity(cells(n),iteration) = sum(AuxImages{1}(labels.Cell==cells(n)));
    CellMeasurements.MedianIntensity(cells(n),iteration) = median(AuxImages{1}(labels.Cell==cells(n)));
end

ModuleDataOut = ModuleData;