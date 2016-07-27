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
    CellMeasurements.MeanIntensity =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.IntegratedIntensity =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MedianIntensity = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.IntensityPercentiles = nan(parameters.TotalCells,parameters.TotalImages,100);
end

% Cycle through each cell and assign measurements
cells = unique(labels.Cell(labels.Cell>0));
for n = 1:length(cells)
    CellMeasurements.MeanIntensity(cells(n),iteration) = mean(AuxImages{1}(labels.Cell==cells(n)));
    CellMeasurements.IntegratedIntensity(cells(n),iteration) = sum(AuxImages{1}(labels.Cell==cells(n)));
    CellMeasurements.MedianIntensity(cells(n),iteration) = median(AuxImages{1}(labels.Cell==cells(n)));
    CellMeasurements.IntensityPercentiles(cells(n),iteration,:) = prctile(AuxImages{1}(labels.Cell==cells(n)),1:100);
end


% Measure cells in 2nd auxiliary image, if it is specified
if ~isempty(AuxImages{2})
    % Mode-balance
    if ~isfield(ModuleData,'distr2')
        [~, ModuleData.distr2] = modebalance(AuxImages{2},2,ModuleData.BitDepth,'measure'); 
    else
        AuxImages{2} = modebalance(AuxImages{2},2,ModuleData.BitDepth,'correct',ModuleData.distr2);
    end


    % On first call, initialize all new CellMeasurements fields 
    if ~isfield(CellMeasurements,'MeanIntensityNuc2')
        % Intensity-based measurement initialization
        CellMeasurements.MeanIntensity2 =  nan(parameters.TotalCells,parameters.TotalImages);
        CellMeasurements.IntegratedIntensity2 =  nan(parameters.TotalCells,parameters.TotalImages);
        CellMeasurements.MedianIntensity2 = nan(parameters.TotalCells,parameters.TotalImages);
        CellMeasurements.IntensityPercentiles2 = nan(parameters.TotalCells,parameters.TotalImages,100);

    end

    % Cycle through each cell and assign measurements
    cells = unique(labels.Cell(labels.Cell>0));
    for n = 1:length(cells)
        CellMeasurements.MeanIntensity2(cells(n),iteration) = mean(AuxImages{2}(labels.Cell==cells(n)));
        CellMeasurements.IntegratedIntensity2(cells(n),iteration) = sum(AuxImages{2}(labels.Cell==cells(n)));
        CellMeasurements.MedianIntensity2(cells(n),iteration) = median(AuxImages{2}(labels.Cell==cells(n)));
        CellMeasurements.IntensityPercentiles2(cells(n),iteration,:) = prctile(AuxImages{2}(labels.Cell==cells(n)),1:100);

    end

end

ModuleDataOut = ModuleData;