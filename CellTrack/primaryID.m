function [output, diagnos] =  primaryID(image0,p, ~)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [output, diagnos] =  primaryID(image0,p,data,~) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PRIMARYID  Identifies candidate regions, in absence of cytoplasmic image, where nuclei can be identified
%
% INPUTS:
% image0        input fluorescent (nuclear) image
% p             parameters struture
%
% OUTPUTS:
% output        all masks needed for tracking purposes
% diagnos       structure with all masks and label matricies (for diagnosis in testing mode)
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%- - - - - - - - - - - - - - - - - - - MAKE "CELL" MASK - - - - - - - - - - - - - - - - - - - - - - -
horizontalEdge = imfilter(image0,fspecial('sobel') /8,'symmetric');
verticalEdge = imfilter(image0,fspecial('sobel')'/8,'symmetric');
nuc_edge = sqrt(horizontalEdge.^2 + verticalEdge.^2);
output.mask1 = nuc_edge>tsaithresh(nuc_edge,false(size(nuc_edge)),2^10);
output.mask1 = bwareaopen(output.mask1,p.NoiseSize);
output.mask_cell =  imdilate(output.mask1,ones(64));

diagnos = output;


