function [CellMeasurements, ModuleData] = fretModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = fretModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% FRETMODULE performs ratiometric FRET measurement in cell nuclei/cytoplasmic regions.
% (FRET image just be in slot #1, CFP image should be in slot #2)
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Get indicies and bwconncomp structures for measurement
iteration  = ModuleData.iter;
nuc_cc = label2cc(labels.Nucleus,0);

% (Use annulus if a cell label matrix isn't defined)
if strcmpi(parameters.ImageType,'none')
    tmp_mask = imdilate(labels.Nucleus>0,diskstrel(parameters.MinNucleusRadius*1.5)); % (Expand by 1.5 radii)
    labels.Cell = IdentifySecPropagateSubfunction(double(labels.Nucleus),zeros(size(labels.Nucleus)),tmp_mask,100);
end

cyto_label = labels.Cell;
cyto_label(labels.Nucleus>0) = 0; % cytoplasm only
cyto_cc = label2cc(cyto_label,0);
cell_cc = label2cc(labels.Cell,0);



% Correct and get ratiometric measurement for the FRET/CFP image pair
try
    fret = AuxImages{1};
    fret = fret - double(parameters.Flatfield{end});
    fret = flatfieldcorrect(fret,double(parameters.Flatfield{1}));
    fret = fret-prctile(fret(:),2); % Background subtract
    fret(fret<16) = 1.6; % add floor to image 

    cfp = AuxImages{2};
    cfp = cfp - double(parameters.Flatfield{end});
    cfp = flatfieldcorrect(cfp,double(parameters.Flatfield{1}));
    cfp = cfp-prctile(cfp(:),2); % Background subtract

    cfp(cfp<16) = 16; % add floor to image 
    fret_image = (fret)./(cfp);
catch me
    % Skip measurement if FRET images are invalid/not found - leave inputs intact
    return;
end

% - - - - NUCLEAR measurements - - - -
% A) Initialize fields
if ~isfield(CellMeasurements,'MeanFRET_nuc')
    CellMeasurements.MeanFRET_nuc =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.IntegratedFRET_nuc =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MedianFRET_nuc = nan(parameters.TotalCells,parameters.TotalImages);

end
% B) Assign measurements
for n = 1:nuc_cc.NumObjects
    CellMeasurements.MeanFRET_nuc(n,iteration) = nanmean(fret_image(nuc_cc.PixelIdxList{n}));
    CellMeasurements.IntegratedFRET_nuc(n,iteration) = nansum(fret_image(nuc_cc.PixelIdxList{n}));
    CellMeasurements.MedianFRET_nuc(n,iteration) = nanmedian(fret_image(nuc_cc.PixelIdxList{n}));
end

% - - - - CYTOPLASMIC/WHOLE-CELL measurements - - - -\
% A) Initialize fields
if ~isfield(CellMeasurements,'MeanFRET_cyto')
    CellMeasurements.MeanFRET_cyto =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.IntegratedFRET_cyto =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MedianFRET_cyto = nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MeanFRET_cell =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.IntegratedFRET_cell =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MedianFRET_cell = nan(parameters.TotalCells,parameters.TotalImages);

    % Grab a "high" measurement - 90th percentile to acct for potential cell shrinkage
    CellMeasurements.HighFRET_cell =  nan(parameters.TotalCells,parameters.TotalImages,10);
    CellMeasurements.HighFRET_cyto =  nan(parameters.TotalCells,parameters.TotalImages,10);

end

% B) Assign measurements
for n = 1:cyto_cc.NumObjects
    CellMeasurements.MeanFRET_cyto(n,iteration) = nanmean(fret_image(cyto_cc.PixelIdxList{n}));
    CellMeasurements.IntegratedFRET_cyto(n,iteration) = nansum(fret_image(cyto_cc.PixelIdxList{n}));
    CellMeasurements.MedianFRET_cyto(n,iteration) = nanmedian(fret_image(cyto_cc.PixelIdxList{n}));

    CellMeasurements.MeanFRET_cell(n,iteration) = nanmean(fret_image(cell_cc.PixelIdxList{n}));
    CellMeasurements.IntegratedFRET_cell(n,iteration) = nansum(fret_image(cell_cc.PixelIdxList{n}));
    CellMeasurements.MedianFRET_cell(n,iteration) = nanmedian(fret_image(cell_cc.PixelIdxList{n}));

    CellMeasurements.HighFRET_cell(n,iteration,:) =  prctile(fret_image(cell_cc.PixelIdxList{n}),linspace(80, 99, 10));
    CellMeasurements.HighFRET_cyto(n,iteration,:) =  prctile(fret_image(cyto_cc.PixelIdxList{n}),linspace(80, 99, 10));
end

