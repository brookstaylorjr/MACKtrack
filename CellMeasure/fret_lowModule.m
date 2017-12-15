function [CellMeasurements, ModuleData] = fret_lowModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = fretModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% FRET_LOWMODULE performs ratiometric FRET measurement in cell nuclei/cytoplasmic regions.
% (FRET image just be in slot #1, CFP image should be in slot #2)
% Flatfields #1 and #2 will be used to correct the FRET/CFP image, respectively.
%
% [Function is identical to fretModule, but uses altered flatfield correction & thresholding options that are 
% better-suited for low-signal images]
%
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

% Parameters for correction/normalization
max_jump = 3;
floor_signal = 32;

% Correct and get ratiometric measurement for the FRET/CFP image pair
try
  
    fret = AuxImages{1};
    fret = fret - double(parameters.Flatfield{end});
    fret = flatfieldcorrect(fret,double(parameters.Flatfield{1}),'subtract');
    fret = fret-prctile(fret(:),2); % Background subtract
    
    cfp = AuxImages{2};
    cfp = cfp - double(parameters.Flatfield{end});
    cfp = flatfieldcorrect(cfp,double(parameters.Flatfield{2}),'subtract');
    cfp = cfp-prctile(cfp(:),2); % Background subtract
    
    % Ensure FRET and CFP images are aligned before point-by-point division
    n = [5 5]; % Total number of image blocks
    [~,r_jumps, c_jumps, maxes] = calculatejump(fret,cfp,n);
    r_jumps(abs(r_jumps)>max_jump) = 0;
    c_jumps(abs(c_jumps)>max_jump) = 0;
    weights = maxes.^2; % (Use exponent to force low-quality correlations to lower weight)
    
    
    % Perform weighted smoothing for each block, then rescale to full image size
    h = fspecial('average');
    c_scaled = imfilter(c_jumps.*weights,h)./imfilter(weights,h);
    r_scaled = imfilter(r_jumps.*weights,h)./imfilter(weights,h);
    c_scaled = round(imresize(c_scaled,size(fret)));
    r_scaled = round(imresize(r_scaled,size(fret)));
    % Apply calculated transformation
    [X, Y] = meshgrid(1:size(fret,2),1:size(fret,1));
    tmp_image = fret;
    pad = max([abs(c_scaled(:)); abs(r_scaled(:))])+1;
    tmp_image = [zeros(pad,size(tmp_image,2)); tmp_image; zeros(pad,size(tmp_image,2))];
    tmp_image = [zeros(size(tmp_image,1),pad), tmp_image, zeros(size(tmp_image,1),pad)];
    fret_remap = tmp_image(sub2ind(size(tmp_image),Y-r_scaled+pad, X - c_scaled+pad));  
    
    fret_image = (fret_remap)./(cfp);    
    fret_image(cfp<floor_signal) = nan; % Omit pixels beneath signal floor
    fret_image(fret_image>10 | fret_image<0) = nan; % Omit spuriously high/low pixels  
    
catch me
    % Skip measurement if FRET images are invalid/not found - leave inputs intact
    return;
end

%%


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
    CellMeasurements.MeanFRET_cell =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.IntegratedFRET_cell =  nan(parameters.TotalCells,parameters.TotalImages);
    
    % Get percentile measurements to see if we can acct for potential cell shrinkage
    CellMeasurements.FRETpctile_cell =  nan(parameters.TotalCells,parameters.TotalImages,15);
    CellMeasurements.FRETpctile_cyto =  nan(parameters.TotalCells,parameters.TotalImages,15);

end

% B) Assign measurements
for n = 1:cyto_cc.NumObjects
    CellMeasurements.MeanFRET_cyto(n,iteration) = nanmean(fret_image(cyto_cc.PixelIdxList{n}));
    CellMeasurements.IntegratedFRET_cyto(n,iteration) = nansum(fret_image(cyto_cc.PixelIdxList{n}));

    CellMeasurements.MeanFRET_cell(n,iteration) = nanmean(fret_image(cell_cc.PixelIdxList{n}));
    CellMeasurements.IntegratedFRET_cell(n,iteration) = nansum(fret_image(cell_cc.PixelIdxList{n}));

    CellMeasurements.FRETpctile_cell(n,iteration,:) =  prctile(fret_image(cell_cc.PixelIdxList{n}),25:5:95);
    CellMeasurements.FRETpctile_cyto(n,iteration,:) =  prctile(fret_image(cyto_cc.PixelIdxList{n}),25:5:95);
end

