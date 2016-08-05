function [CellMeasurements, ModuleData] = tmrmModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% TMRMMODULE measures cellular intensity in AuxImages{1} channel
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
if ~isfield(ModuleData,'distr2')
    [AuxImages{1}, ModuleData.distr2] = modebalance(AuxImages{1},2,ModuleData.BitDepth,'measure');
    ModuleData.tmrm_thresh = quickthresh(AuxImages{1},false(size(AuxImages{1})),'log');

else
    AuxImages{1} = modebalance(AuxImages{1},2,ModuleData.BitDepth,'correct',ModuleData.distr2);
end

    
tmrm = AuxImages{1};
mask1 = tmrm > ModuleData.tmrm_thresh;

% Initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'TMRM_cytoMean')
    % Intensity-based measurement initialization
    CellMeasurements.TMRM_cytoMean =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.TMRM_cytoMedian =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.TMRM_nucMean =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.TMRM_nucMedian =  nan(parameters.TotalCells,parameters.TotalImages);    
end

% Cycle through each cell and assign measurements
cells = unique(labels.Cell(labels.Cell>0));
for n = 1:length(cells)
    nucleus = labels.Nucleus==cells(n);
    cytoplasm = (labels.Cell==cells(n)) &~nucleus;
    
    CellMeasurements.TMRM_cytoMean(cells(n),iteration) = mean(tmrm(mask1 & cytoplasm));
    CellMeasurements.TMRM_cytoMedian(cells(n),iteration) = median(tmrm(mask1 & cytoplasm));
    CellMeasurements.TMRM_nucMean(cells(n),iteration) = mean(tmrm(mask1 & nucleus));
    CellMeasurements.TMRM_nucMedian(cells(n),iteration) = mean(tmrm(mask1 & nucleus));
end

ModuleDataOut = ModuleData;