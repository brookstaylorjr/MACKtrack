function [CellMeasurements, ModuleData] = fret_nucModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellMeasurements, ModuleData] = fret_nucModule(CellMeasurements, parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% FRET_NUCMODULE performs ratiometric FRET measurement in cell nuclei/cytoplasmic regions.
% (FRET image just be in slot #1, CFP image should be in slot #2).
% Nuclear-only signal is assumed, so a more aggressive normalization schema is used.
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


% Define background region
drop_region = imdilate(labels.Nucleus>0,diskstrel(2*parameters.MinNucleusRadius));

% Correct and get ratiometric measurement for the FRET/CFP image pair
try
    fret = AuxImages{1};
    fret = fret - double(parameters.Flatfield{end});
    fret = flatfieldcorrect(fret,double(parameters.Flatfield{1}),'subtract');
    % Fit unimodal normal distribution to image background. Subtract that mean from original image
    [~, d] = modebalance(fret(~drop_region),1,ModuleData.BitDepth,'measure');
    fret = fret-d(1);
    fret(fret<16) = 1.6; % add floor to image 

    
    cfp = AuxImages{2};
    cfp = cfp - double(parameters.Flatfield{end});
    cfp = flatfieldcorrect(cfp,double(parameters.Flatfield{1}),'subtract');
    % Fit unimodal normal distribution to image background. Subtract that mean from original image
    [~, d] = modebalance(cfp(~drop_region),1,ModuleData.BitDepth,'measure');
    cfp = cfp-d(1);
    cfp(cfp<16) = 16; % add floor to image
    
    
    % Ensure FRET and CFP images are aligned - use normal image alignment tool
    n = [5 5];
    [~,r_jumps, c_jumps, maxes] = calculatejump(fret,cfp,n);
    weights = maxes.^2; % Rescale maxes (force low-quality correlations to near-zero weight)
    h = fspecial('average',[3 3]);
    c_scaled = imfilter(c_jumps.*weights,h)./imfilter(weights,h);
    r_scaled = imfilter(r_jumps.*weights,h)./imfilter(weights,h);
    c_scaled = round(imresize(c_scaled,size(fret)));
    r_scaled = round(imresize(r_scaled,size(fret)));
    [X, Y] = meshgrid(1:size(fret,2),1:size(fret,1));
    tmp_image = cfp;
    pad = max([abs(c_scaled(:)); abs(r_scaled(:))])+1;
    tmp_image = [zeros(pad,size(tmp_image,2)); tmp_image; zeros(pad,size(tmp_image,2))];
    tmp_image = [zeros(size(tmp_image,1),pad), tmp_image, zeros(size(tmp_image,1),pad)];
    idx = sub2ind(size(tmp_image),Y+r_scaled+pad, X + c_scaled+pad);
    a = tmp_image(idx);
    
    fret_image = (fret)./(a);
catch me
    % Skip measurement if FRET images are invalid/not found - leave inputs intact
    return;
end

% - - - - NUCLEAR measurements - - - -
% A) Initialize fields
if ~isfield(CellMeasurements,'MeanFRET_nuc')
    CellMeasurements.MeanFRET_nucAlt =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.IntegratedFRET_nucAlt =  nan(parameters.TotalCells,parameters.TotalImages);
    CellMeasurements.MedianFRET_nucAlt = nan(parameters.TotalCells,parameters.TotalImages);

end
% B) Assign measurements
for n = 1:nuc_cc.NumObjects
    CellMeasurements.MeanFRET_nucAlt(n,iteration) = nanmean(fret_image(nuc_cc.PixelIdxList{n}));
    CellMeasurements.IntegratedFRET_nucAlt(n,iteration) = nansum(fret_image(nuc_cc.PixelIdxList{n}));
    CellMeasurements.MedianFRET_nucAlt(n,iteration) = nanmedian(fret_image(nuc_cc.PixelIdxList{n}));
end
