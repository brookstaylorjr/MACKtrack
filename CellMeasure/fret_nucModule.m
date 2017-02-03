function [CellMeasurements, ModuleData] = fret_nucModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = fret)nucModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% FRET_NUCMODULE 
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

iteration  = ModuleData.iter;
tmp_label = labels.Nucleus;
measure_cc = label2cc(tmp_label,0);


% Get ratiometric measurement for FRET
fret = AuxImages{1}; 
fret = flatfieldcorrect(fret,double(parameters.Flatfield{1}));
fret = fret- prctile(fret(:),2);
fret(fret<20) = 20;
cfp = AuxImages{2};
cfp = flatfieldcorrect(cfp,double(parameters.Flatfield{1}));
cfp = cfp - prctile(cfp(:),2);
cfp(cfp<20) = 20;
fret_image = (fret)./(cfp);


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
