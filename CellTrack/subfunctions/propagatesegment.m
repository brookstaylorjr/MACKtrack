function [label_out] = propagatesegment(seeds, mask, image_in, strel_size, nuclei, lambda)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PROPAGATESEGEMENT performs 'propagate' algorithm (taken from CellProfiler), then corrects
% small errors
%
% seeds        label matrix of seed objects
% mask         mask of areas to fill/segment 
% image_in     original image
% strel_size   size of structuring element used to identify thin pieces of cells
% nuclei       nuclear seeds (default to "seeds" if not specified)
% lambda       (optional) scaling parameter for Propagate. Default is 0.02
%
% label_out    final image
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin<6
    lambda = 0.02;
    if nargin<5
        nuclei = seeds;
    end
end

image_clamp = abs((image_in-prctile(image_in(:),0.02))/diff(prctile(image_in(:),[0.02 98])));
image_clamp(image_clamp<0) = 0; image_clamp(image_clamp>1) = 1;

% Run propagate (requires compiled binary)
label1 = IdentifySecPropagateSubfunction(double(seeds),double(image_clamp),mask,lambda);

%%
% Correct for errors segmentation: look for thin pieces
SEround = diskstrel(strel_size);
borders = (imdilate(label1,ones(3))-label1)>0.5;
thin_pieces = label1>0;
thin_pieces(borders) = 0;
thin_pieces = (thin_pieces)&~imopen(thin_pieces, SEround);

% Remove small bits from opened image
thin_pieces = bwareaopen(thin_pieces,4);
% Add back border that was subtracted earlier
thin_pieces = thin_pieces|(imdilate(thin_pieces,ones(3)) & borders & (label1>0) );

% Now, reassign objects based on who they share longest border with.
label_dilated = label1;
label_dilated(thin_pieces) = 0;
label_dilated = imdilate(label_dilated,ones(3));
thin_cc = bwconncomp(thin_pieces);
seed_reassign = label1;

for i = 1:thin_cc.NumObjects
    newObj = label_dilated(thin_cc.PixelIdxList{i});
    newObj(newObj==0) = [];
    newObj = mode(newObj);
    seed_reassign(thin_cc.PixelIdxList{i}) = newObj;
end
seed_reassign(isnan(seed_reassign)) = 0;

% Make sure no "islands" exist- all cells should be contiguous with original nuclei 
seed_contig = seed_reassign;
erode1 = imerode(seed_contig,ones(3));
dilate1 = imdilate(seed_contig,ones(3));
dilate1(seed_contig==0) = 0;
borders1 = seed_contig - erode1;
borders2 = dilate1 - seed_contig;
seed_contig( (borders1~=0) | (borders2~=0) )=0;
seed_contig = bwlabel(seed_contig>0,4);
seed_contig = removemarked(seed_contig,nuclei>0,'keep');
label2 = IdentifySecPropagateSubfunction(double(nuclei),double(image_clamp),seed_contig>0,lambda);


% Do final propagation
label_out = IdentifySecPropagateSubfunction(double(label2),double(image_clamp),mask,lambda);
label_out(isnan(label_out)) = 0;
