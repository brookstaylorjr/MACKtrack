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
    thresh = 16;
    fret  = AuxImages{1};   
    fret = flatfieldcorrect(fret-double(ff{end}),double(ff{1}));
    fret = fret- prctile(fret(:),2);
    fret(fret<thresh) = thresh/10;

    cfp  = AuxImages{2};   
    cfp = flatfieldcorrect(cfp-double(ff{end}),double(ff{1}));
    cfp = cfp - prctile(cfp(:),2);
    cfp(cfp<thresh) = thresh;
    
    % Ensure FRET and CFP images are aligned before point-by-point division
    max_jump = 3;
    n = [5 5];
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
catch me
    % Skip measurement if FRET images are invalid/not found - leave inputs intact
    return;
end

% - - - - NUCLEAR measurements - - - -
% A) Initialize fields
if ~isfield(CellMeasurements,'MeanFRET_nucAlt')
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
