function filter_out = gauss2D(stdev_in)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% GAUSS    is a simple function which computes a Gaussian filter based on the input
% standard deviation
%
% stdev_in     input standard deviation
%
% filter_out   output gaussian filter
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
filter_out = fspecial('gaussian',6*ceil(stdev_in)+1,stdev_in);