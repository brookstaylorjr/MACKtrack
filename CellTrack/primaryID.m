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

% Shrink image if it's very large (>2 megapixels, e.g. from nonbinned image)
orig_size = size(image0);
scale_factor = 1000/orig_size(1);
n = round(log(1/scale_factor)/log(2));
scale_factor = 1/(2^n); % Round to nearest pwr of 2
if numel(image0) > 2e6
    image0 = imresize(image0,scale_factor,'nearest');
end
nuc_smooth = imfilter(image0,gauss2D(scale_factor*p.MinNucleusRadius/2),'replicate');
stdev_in = round(scale_factor*p.MinNucleusRadius);
h = fspecial('log',8*ceil(stdev_in)+1, stdev_in);
nuc_edge = abs(imfilter(nuc_smooth, h,'symmetric'));

output.mask1 = nuc_edge>exp(trithreshold(log(nuc_edge)));
output.mask1 = imresize(output.mask1,orig_size,'nearest');
output.mask_cell =  imdilate(output.mask1,ones(ceil(p.MinNucleusRadius)));

diagnos = output;
