function [CellMeasurements, ModuleData] = fretModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = fretModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% FRETMODULE 
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

iteration  = ModuleData.iter;
tmp_label = labels.Cell;
tmp_label(labels.Nucleus>0) = 0;
measure_cc = label2cc(tmp_label,0);

[~,dist1] = modebalance(AuxImages{1},2,16,'measure');
AuxImages{1} = AuxImages{1}-dist1(1);

[~,dist2] = modebalance(AuxImages{2},2,16,'measure');
AuxImages{2} = AuxImages{2}-dist2(1);

AuxImages{1} = AuxImages{1}-(min(AuxImages{1}(:)))+5;
AuxImages{2} = AuxImages{2} - min(AuxImages{2}(:))+5;
fret_image = (AuxImages{1})./(AuxImages{2});


% On first call, initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'MeanFRET')
    % Intensity-based measurement initialization
    CellMeasurements.MeanFRET =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.IntegratedFRET =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MedianFRET = nan(parameters.TotalCells,parameters.TotalImages);
end

% Cycle through each cell and assign measurements
for n = 1:measure_cc.NumObjects
    CellMeasurements.MeanFRET(n,iteration) = nanmean(fret_image(measure_cc.PixelIdxList{n}));
    CellMeasurements.IntegratedFRET(n,iteration) = nansum(fret_image(measure_cc.PixelIdxList{n}));
    CellMeasurements.MedianFRET(n,iteration) = nanmedian(fret_image(measure_cc.PixelIdxList{n}));
end
