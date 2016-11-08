function [output, diagnos] =  fluorescenceID(cell_image,p, X)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [output, diagnos] =  fluorescenceID(image0,p,data,~) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% FLUORESCENCEID  Identifies candidate regions for cells (foreground) in a fluorescent image
%
% INPUTS:
% image0        input fluorescent (cell) image
% p             parameters struture
% image1        input fluorescent (nucleus) image

%
% OUTPUTS:
% output        all masks needed for tracking purposes
% diagnos       structure with all masks and label matricies (for diagnosis in testing mode)
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%- - - - - - - - - - - - - - - - - - - MAKE "CELL" MASK - - - - - - - - - - - - - - - - - - - - - - -
%%
warning off MATLAB:nearlySingularMatrix
pStar = (X'*X)\(X')*double(cell_image(:));
% Apply correction
cell_corr = reshape((double(cell_image(:) - X*pStar)),size(cell_image));
cell_corr = cell_corr-min(cell_corr(:)); % Set minimum to zero

%%
[~,bg_dist] = modebalance(cell_image, 1, 16, 'measure');


%nuc_filt = imfilter(image0,gauss2D(p.MinNucleusRadius/2),'replicate'); % Gaussian filtered
% Combine two threshold variants (Otsu threshold + MoG threshold)
thresh1 = min([ quickthresh(cell_image,false(size(cell_image)),'none'), bg_dist(1)+2*bg_dist(2)]);
diagnos.search = cell_image>thresh1;




%%

diagnos.search_dilate = imdilate(diagnos.search,ones(floor(p.MaxNucleusRadius*1.5)));
output.mask1 = diagnos.search_dilate;
output.mask_cell = diagnos.search_dilate;



