function [CellMeasurements, ModuleDataOut] = spotcountModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SPOTCOUNTMODULE applies a high-pass filter to image, then counts regional maxima
% It requires a 2ns Aux image to be passed- this should be a in image with
% lots of spots, you build a threshold off of.
%
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
iteration = ModuleData.iter;
spotSize = 5;

% On 1st iteration, threshold using the secondarily-passed image
if iteration ==1
    spots_hp = AuxImages{2} - imfilter(AuxImages{2},gauss2D(spotSize/2),'replicate');
    % Threshold- use noisethreshold method on image intensity, using a v. small noise size (2)
    search_range = prctile(spots_hp(:),[75 99.99]); % Assumes spot image is pretty sparse
    ModuleData.spotthresh = noisethresh(spots_hp,false(size(spots_hp)),search_range,3);
    % Initialize data
    CellMeasurements.NumSpots = zeros(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.Spots = cell(parameters.TotalCells,parameters.TotalImages);
end

spots_hp = AuxImages{1} - imfilter(AuxImages{1},gauss2D(spotSize/2),'replicate');
mask1 = spots_hp > ModuleData.spotthresh;
% Filter for spots that have local maxima centered near center of spot
test = imfilter(AuxImages{1},gauss2D(1),'replicate');
mask2 = imdilate((imregionalmax(test)&mask1),ones(5)) & mask1;
% Filter speckle noise and label
spot_label = bwlabel(bwareaopen(mask2,3));
cells = unique(labels.Cell); cells(cells==0) = [];
cell_cc = label2cc(labels.Cell);
spot_measurements = cell2mat(struct2cell(regionprops(spot_label,AuxImages{1},'MeanIntensity')));

for n =1:length(cells)
    spot_ids = unique(spot_label(cell_cc.PixelIdxList{n}));spot_ids(spot_ids==0) = [];
    CellMeasurements.Spots{cells(n),iteration} = spot_measurements(spot_ids);
    CellMeasurements.NumSpots(cells(n),iteration) = length(spot_ids);
end

ModuleDataOut = ModuleData;