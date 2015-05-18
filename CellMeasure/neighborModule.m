function [CellMeasurements, ModuleDataOut] = neighborModule(CellMeasurements,parameters, labels, ~, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%  [CellMeasurements, ModuleDataOut] = neighborModule(CellMeasurements,parameters, labels, ~, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NEIGHBORMODULE  measure how many cells each cell is touching, as well as how many are nearby
%
% INPUT:
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement
%
% OUTPUT:
% CellMeasurements    structure with fields corresponding to cell measurements
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Default conversion factor for my images = 4.03 pix/micron
conversion = 4.03;



% On first call, initialize all new CellMeasurements fields
if ~isfield(CellMeasurements,'OccupiedPerimeter')
    CellMeasurements.OccupiedPerimeter = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.Neighbor10 =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.Neighbor25 =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.Neighbor50 =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.Neighbor100 =  nan(parameters.TotalCells,parameters.TotalImages);

end

% Pull cells from label matrix, convert to conncomp, and get centroids
cells = unique(labels.Cell(labels.Cell>0));
cc1 = label2cc(labels.Cell);
props = regionprops(cc1, 'Centroid');
tmpcell = struct2cell(props);
tmpmat = cell2mat(tmpcell(1,:));
centroidx = tmpmat(1:2:end);
centroidy = tmpmat(2:2:end);

% Occupied Perimeter: percentage of cell perimeter that touches other cells
% Start by getting each cell's perimeter pixels
erode1 = imerode(labels.Cell,ones(3));
dilate1 = imdilate(labels.Cell,ones(3));
dilate1(labels.Cell==0) = 0;
borders1 = labels.Cell - erode1;
borders2 = dilate1 - labels.Cell;
border_mask =  (borders1~=0) | (borders2~=0);
% Drop perimeter pixels that border background
touch_mask = border_mask;
touch_mask(erode1==0) = 0;
% Measure the ratio between the two
for n = 1:length(cells)
    CellMeasurements.OccupiedPerimeter(cells(n),ModuleData.iter) = sum(touch_mask(cc1.PixelIdxList{n})) / sum(border_mask(cc1.PixelIdxList{n}));
end

% Nearby cells: count number of cells within a certain cistance of starting cells
% Pair cells, and look at (discretized) distance between their nearest edges
all_dist = nan(length(cells),length(cells));
for m = 1:(length(centroidx)-1)
    for n = (m+1):length(centroidx)
        % Calculate diagonal line connecting cell centroids
        dist1 = max(abs(centroidx(n)-centroidx(m)) , abs(centroidy(n)-centroidy(m)))+2;
        r = round(linspace(centroidy(n),centroidy(m),dist1));
        c = round(linspace(centroidx(n),centroidx(m),dist1));
        line_ind = sub2ind(cc1.ImageSize,r,c);
        % Remove any line_ind that is member of either cell in pair
        line_ind(ismember(line_ind,[cc1.PixelIdxList{n};cc1.PixelIdxList{m}])) = [];
        % Assign resultant line length into results matrix
        all_dist(m,n) = length(line_ind);
        all_dist(n,m) = length(line_ind);
    end   
end

% Count cells beneath distance thresholds - use 4.03px/um conversion factor
for n=1:length(cells)
    % Get distance from cell to image edge- don't assign any measurements greater than this value
    [r, c] = ind2sub(size(labels.Cell), cc1.PixelIdxList{n});
    edgedist = min([size(labels.Cell,1)-max(r), size(labels.Cell,2)-max(c), min(r)-1,min(c)-1]);
    if edgedist > (10*conversion)
        CellMeasurements.Neighbor10(cells(n),ModuleData.iter) =  sum(all_dist(:,n) < (10*conversion));
    end
    if edgedist > (25*conversion)
        CellMeasurements.Neighbor25(cells(n),ModuleData.iter) =  sum(all_dist(:,n) < (25*conversion));
    end
    if edgedist > (50*conversion)
        CellMeasurements.Neighbor50(cells(n),ModuleData.iter) =  sum(all_dist(:,n) < (50*conversion));
    end
    if edgedist > (100*conversion)
        CellMeasurements.Neighbor100(cells(n),ModuleData.iter) =  sum(all_dist(:,n) < (100*conversion));
    end
end

ModuleDataOut= ModuleData;
