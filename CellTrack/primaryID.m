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

%[~,bg_dist] = modebalance(image0, 2, 16, 'measure');


%nuc_filt = imfilter(image0,gauss2D(p.MinNucleusRadius/2),'replicate'); % Gaussian filtered
% Combine two threshold variants (Otsu threshold + MoG threshold)
%thresh1 = min([quickthresh(image0,false(size(image0)),'none'), bg_dist(1)+4*bg_dist(2)]);



%diagnos.search = image0>thresh1;
%diagnos.search_dilate = imdilate(diagnos.search,ones(floor(p.MaxNucleusRadius*1.5)));
%output.mask1 = diagnos.search_dilate;
%output.mask_cell = diagnos.search_dilate;

output.mask1 = true(size(image0));
output.mask_cell =  true(size(image0));

diagnos = output;


