function [CellMeasurements, ModuleDataOut] = edgeModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% EDGEMODULE  measures edge content within cell mask (designed to be used on DIC channel)
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current ModuleData.iter, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


% On first call, initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'MedianEdge')
    % Intensity-based measurement initialization
    CellMeasurements.MedianEdge = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MedianEdge_erode3 = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MedianEdge_erode5 = nan(parameters.TotalCells,parameters.TotalImages);
end

% Get cells in each frame
cells = unique(labels.Cell);
cells(cells==0) = [];

% Convert image to edge magnitude
horizontalEdge = imfilter(AuxImages{1},fspecial('sobel') /8,'symmetric');
verticalEdge = imfilter(AuxImages{1},fspecial('sobel')'/8,'symmetric');
EdgeMag = sqrt(horizontalEdge.^2 + verticalEdge.^2);


% Cycle through each image - get mask of cell, measure median edge content and assign measurements
for i = 1:length(cells)
    cell_mask = labels.Cell==cells(i);
    % Assign measurements
    CellMeasurements.MedianEdge(cells(i),ModuleData.iter) = median(EdgeMag(cell_mask));
    CellMeasurements.MedianEdge_erode3(cells(i),ModuleData.iter) = median(EdgeMag(imerode(cell_mask),ones(3)));
    CellMeasurements.MedianEdge_erode5(cells(i),ModuleData.iter) = median(EdgeMag(imerode(cell_mask),ones(5)));
end

ModuleDataOut = ModuleData;
