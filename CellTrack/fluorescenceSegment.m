function [output, diagnos] = fluorescenceSegment(data, image_in, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [output, diagnos] = fluorescenceSegment(data, original, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% FLUORESECNCESEGMENT 
%
%
% INPUT:
% data       tracking info structure: cell mask, nucleus label matrix, and supporting information
% image_in    original cell image
% p           parameters structure from SetupTracking.m
%
% OUTPUT:
% output      all information (masks) needed for tracking
% diagnos     major masks/thresholds created as intermediates
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

diagnos = struct;
output.img_straight = image_in;
output.cells = propagatesegment(data.nuclei, data.mask_cell>0, image_in,...
    round(p.MinCellWidth/2),data.nuclei,0.02);

% Fill all holes that are internal to a single object
cell_holes = imdilate(imfill(output.cells>0,'holes') &~ (output.cells),ones(3));
hole_cc = bwconncomp(cell_holes);
for i = 1:hole_cc.NumObjects
    hole_vals = unique(output.cells(hole_cc.PixelIdxList{i}));
    hole_vals(hole_vals==0) = [];
    if numel(hole_vals) == 1
        output.cells(hole_cc.PixelIdxList{i}) = hole_vals;
    end
end


% Save all information under diagnostic struct
diagnos = combinestructures(diagnos,output);