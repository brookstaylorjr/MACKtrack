function [label_out, data] = phaseSegment(data, image_in, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% phaseSegment   Segment clusters of cells in the binary mask of a phase contrast image.
% (Based on CellProfiler's "propogate" algorithm, but improves result using contour points
% and intensity mask.)
%
%
% queue       tracking info structure: cell mask, nucleus label matrix, and supporting information
% image_in    image_in phase contrast image
% p           parameters structure from SetupTracking.m
%
% label_out   label matrix showing one cell per nucleus in data.nuclei
% data        diagnostic information
%
% Subfunctions
% matchclosest.m, IdentifySecPropagateSubfunction.cpp (compiled w/ mex)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%% Setup: create basic masks/label matrices from phase/nuclear images
cell_mask = data.mask_cell>0;
nuc_mask = data.nuclei>0;
halos = data.halos;
SEround = strel('disk',floor(p.MinCellWidth/2),4);

%% Inflection point splitting

% Find all candidate inflection points
diagnos.modifier_mask = findinflection(cell_mask,nuc_mask,p);
if nargout>1
    diagnos.all_inflection = diagnos.modifier_mask+nuc_mask+cell_mask;
end
%% Do final propagation/segementation, and correction.
% Turn off pixels from modifier mask, do initial segmentation
data.modifier_mask = removemarked(bwlabel(data.modifier_mask),data.nuclei>0)> 0;
data.modifier_mask = imdilate(data.modifier_mask,ones(3));
cell_mask(data.modifier_mask) = 0;

lambda = .02; % (Has very little effect on final result)

if strcmp(p.ImageType,'phase')
    % Use halos mask to possibly identify isolated cells (surrounded by halo)
    darkCells = (~halos)& cell_mask;
    darkCells = imopen(darkCells,SEround);
    darkLabeled = bwconncomp(darkCells);
    data.seed_label1 = data.nuclei;
    for i = 1:darkLabeled.NumObjects
        a = unique(data.nuclei(darkLabeled.PixelIdxList{i}));
        a(a==0) = [];
        if length(a)==1
            data.seed_label1(darkLabeled.PixelIdxList{i}) = a;
        end
    end
    data.seed_label2 = IdentifySecPropagateSubfunction(double(data.seed_label1),double(image_in),cell_mask,lambda);
    data.seed_label2(isnan(data.seed_label2)) = 0;
end

% b) Correct for errors in round-cell segmentation.
borders = (imdilate(data.seed_label2,ones(3))-data.seed_label2)>0.5;
thin_pieces = data.seed_label2>0;
thin_pieces(borders) = 0;
thin_pieces = (thin_pieces)&~imopen(thin_pieces, SEround);
% Remove small noise from opened image
thin_pieces = bwareaopen(thin_pieces,4);
% Add back border that was subtracted earlier
thin_pieces = thin_pieces|(imdilate(thin_pieces,ones(3)) & borders & (data.seed_label2>0) );
% Now, reassign objects based on who they share longest border with.
label_dilated = data.seed_label2;
label_dilated(thin_pieces) = 0;
label_dilated = imdilate(label_dilated,ones(3));
reassignLabel = bwconncomp(thin_pieces);
data.seed_reassign = data.seed_label2;
for i = 1:reassignLabel.NumObjects
    newObj = label_dilated(reassignLabel.PixelIdxList{i});
    newObj(newObj==0) = [];
    newObj = mode(newObj);
    data.seed_reassign(reassignLabel.PixelIdxList{i}) = newObj;
end
data.seed_reassign(isnan(data.seed_reassign)) = 0;

% Make sure no "islands" exist- all cells should be contiguous with their nuclei
data.seed_contig = data.seed_reassign;
data.seed_contig((imdilate(data.seed_contig,ones(5))-data.seed_contig)> 0) = 0; % Turns borders off
data.seed_contig = bwlabel(data.seed_contig>0,4);
data.seed_contig = removemarked(data.seed_contig,data.nuclei>0,'keep');
data.seed_label3 = IdentifySecPropagateSubfunction(double(data.nuclei),double(image_in),data.seed_contig>0,lambda);

% Do final assignment
label_out = IdentifySecPropagateSubfunction(double(data.seed_label3),double(image_in),data.mask_cell,lambda);
label_out(isnan(label_out)) = 0;