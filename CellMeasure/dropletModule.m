function [CellMeasurements, ModuleData] = dropletModule(CellMeasurements,parameters, labels, AuxImages, ModuleData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% DROPLETMODULE computes the amount of lipid droplets (in a BF image) in a given cell
%
% CellMeasurements    structure with fields corresponding to cell measurements
%
% parameters          experiment data (total cells, total images, output directory)
% labels              Cell,Nuclear label matricies (labels.Cell and labels.Nucleus)
% AuxImages           images to measure
% ModuleData          extra information (current iteration, etc.) used in measurement 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% If cells were not segmented; do annulus to estimate cytoplasmic area.
if strcmpi(parameters.ImageType,'none')
    tmp_mask = imdilate(labels.Nucleus>0,diskstrel(parameters.MinNucleusRadius*3)); % (Expand by 3 radii bc FAT CELLS ARE HUUUGE)
    labels.Cell = IdentifySecPropagateSubfunction(double(labels.Nucleus),zeros(size(labels.Nucleus)),tmp_mask,100);
end




% Get indicies and bwconncomp structures for measurement
iteration  = ModuleData.iter;

% Make conncomp structures for nucleus, cell, and cytoplasm
nuc_cc = label2cc(labels.Nucleus,0);
cell_cc = label2cc(labels.Cell,0);
cyto_cc = cell_cc;
for i = 1:cyto_cc.NumObjects
    if ~isempty(cyto_cc.PixelIdxList{i})
        cyto_cc.PixelIdxList{i}(ismember(cyto_cc.PixelIdxList{i},nuc_cc.PixelIdxList{i})) = [];
    end
end



[circ_centers, circ_rad] = imfindcircles(AuxImages{1},[4 15],'Method','TwoStage','Sensitivity',0.9,'EdgeThreshold',0.15);
accum_img = false(size(AuxImages{1}));
for i = 1:length(rad)
    tmp = false(size(AuxImages{1}));
    tmp(round(circ_centers(i,2)),round(circ_centers(i,1))) = 1;
    tmp = imdilate(tmp,diskstrel(circ_rad(i)));
    
    accum_img = accum_img | tmp;
end


% A) Initialize fields
if ~isfield(CellMeasurements,'LipidLike')
    CellMeasurements.LipidLike =  nan(parameters.TotalCells,parameters.TotalImages);
end
    % B) Assign measurements
    for n = 1:cyto_cc.NumObjects
        CellMeasurements.LipidLike(n,iteration) = nansum(accum_img(cyto_cc.PixelIdxList{n}));
    end





