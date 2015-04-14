function [CellMeasurements, ModuleDataOut] = morphologyModule(CellMeasurements,parameters, labels, ~, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% morphologyModule measures simple aspects of nuclei/cellular morphologoy and movement
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
iteration = ModuleData.iter;

% On first call, initialize all new CellMeasurements fields 
if ~isfield(CellMeasurements,'Area')
    % Morphology-based cell measurements
    CellMeasurements.Area = zeros(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.Compactness = zeros(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.AxisRatio = zeros(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.Perimeter = zeros(parameters.TotalCells,parameters.TotalImages);
    % Movement measurements
    CellMeasurements.MovementCell = zeros(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MovementNucleus = zeros(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.CellChange = zeros(parameters.TotalCells,parameters.TotalImages);
    % Track each cell's location
    CellMeasurements.CentroidX = zeros(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.CentroidY = zeros(parameters.TotalCells,parameters.TotalImages);
    % Create tmp storage files for multiple time-point measurements
    CellMeasurements.tmp.CentroidCell = zeros(parameters.TotalCells,2);
    CellMeasurements.tmp.CentroidNuc = zeros(parameters.TotalCells,2);
    CellMeasurements.tmp.PrevLabel = [];
end

% Get unique cells, conncomp structure, and regionprops
cells = unique(labels.Cell(labels.Cell>0));
cc1 = label2cc(labels.Cell);
cellProps = regionprops(labels.Cell,'Area','Perimeter','MajorAxisLength','MinorAxisLength','Centroid');
nucProps = regionprops(labels.Nucleus,'Centroid');


% Cycle through each image and assign measurements
for n = 1:length(cells)
    % Basic morpholody measurements
    CellMeasurements.Area(cells(n),iteration) = cellProps(cells(n)).Area;
    CellMeasurements.Perimeter(cells(n),iteration) = cellProps(cells(n)).Perimeter;
    CellMeasurements.AxisRatio(cells(n),iteration) = cellProps(cells(n)).MajorAxisLength/cellProps(cells(n)).MinorAxisLength;
    CellMeasurements.Compactness(cells(n),iteration) = ((cellProps(cells(n)).Perimeter).^2) / (4*pi*cellProps(cells(n)).Area);
    % Cell movement measurements - only valid after cell has existed for 1 frame
    if max(CellMeasurements.tmp.CentroidCell(cells(n),:)) > 0
        CellMeasurements.MovementCell(cells(n),iteration) = norm(cellProps(cells(n)).Centroid-CellMeasurements.tmp.CentroidCell(cells(n),:));
        CellMeasurements.MovementNucleus(cells(n),iteration) =  norm(nucProps(cells(n)).Centroid-CellMeasurements.tmp.CentroidNuc(cells(n),:));   
        overlap = sum(CellMeasurements.tmp.PrevLabel(cc1.PixelIdxList{n})==cells(n))/length(cc1.PixelIdxList{n});
        CellMeasurements.CellChange(cells(n),iteration) = 1-overlap;
    end
    % Location measurements
    CellMeasurements.CentroidX(cells(n),iteration) = cellProps(cells(n)).Centroid(1);
    CellMeasurements.CentroidY(cells(n),iteration) = cellProps(cells(n)).Centroid(2);
    % Update centroid measurements in memory
    CellMeasurements.tmp.CentroidCell(cells(n),:) = cellProps(cells(n)).Centroid;
    CellMeasurements.tmp.CentroidNuc(cells(n),:) = nucProps(cells(n)).Centroid;
end
% Update label matrix in memory
CellMeasurements.tmp.PrevLabel = labels.Cell;
ModuleDataOut = ModuleData;