function xys_out = xls(xy_array, sites)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% xys_out = xls(xy_array, sites)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% XLS generates XY positions from a specified array of wells (e.g. E04 -> 504), and a number of
% sites per well.
%
% Example: xls([504 505], 2) generates [504.01, 504.02, 505.01, 505.02]
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

xys_out = [];

for i = 1:length(xy_array)
    xys_out = cat(2,xys_out,xy_array(i)+((1:sites)/100));
end