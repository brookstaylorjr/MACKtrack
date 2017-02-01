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
    CellMeasurements.Area = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.Compactness = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.AxisRatio = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.Perimeter = nan(parameters.TotalCells,parameters.TotalImages);
    % Track each cell's location
    CellMeasurements.CentroidX = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.CentroidY = nan(parameters.TotalCells,parameters.TotalImages);
    
    % Create tmp storage files for multiple time-point measurements
    CellMeasurements.tmp.CentroidCell = nan(parameters.TotalCells,2);
    CellMeasurements.tmp.CentroidNuc = nan(parameters.TotalCells,2);
    CellMeasurements.tmp.PrevLabel = [];
end

    
if (ModuleData.iter>1) && ~isfield(CellMeasurements,'MovementCell')
    % Movement measurements
    CellMeasurements.MovementCell = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MovementNucleus = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.CellChange = nan(parameters.TotalCells,parameters.TotalImages);
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
    if (ModuleData.iter>1)
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