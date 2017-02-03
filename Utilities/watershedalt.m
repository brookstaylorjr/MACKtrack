function [labelOut] = watershedalt(imageIn,maskIn, connectivity)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [labelOut] = watershedalt(imageIn,maskIn, connectivity)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% WATERSHEDALT fluorescence image variant on watershed- conditions imageIn using maskIn, 
% then modifies output to ensure borders and background are set to zero.
%
% imageIn         base fluorescence image to be watershed-transformed
% maskin          limiting binary mask, defines foreground/background
% connectivity    specify 4 or 8 connectivity- defaults to 4.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if nargin < 3
    connectivity = 4;
end

% Drop background to -Inf/ invert image
imageIn(~maskIn) = Inf;
labelOut = watershed(-imageIn,connectivity);


% Drop any background-identical objects
for i = unique(labelOut(~maskIn))'
    labelOut(labelOut==i) = 0;
end

% Relabel contiguously
labelOut = bwlabel(labelOut>0,4);