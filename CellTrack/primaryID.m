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
%%

nuc_smooth = imfilter(image0,gauss2D(p.MinNucleusRadius/2),'replicate');
stdev_in = p.MinNucleusRadius;
h = fspecial('log',8*ceil(stdev_in)+1, stdev_in);
nuc_edge = abs(imfilter(nuc_smooth, h,'symmetric'));

output.mask1 = nuc_edge>exp(trithreshold(log(nuc_edge)));
output.mask_cell =  imdilate(output.mask1,ones(ceil(p.MinNucleusRadius)));

diagnos = output;
