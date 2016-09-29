function [output, diagnos] = primarySegment(data, image_in, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [output, diagnos] = primarySegment(data, original, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% PRIMARYSEGMENT is a dummy function - it's akin to dic or phaseSegment, but in this case 
% there are no cells to segment, so "cells" are just copied from nuclei.
%
% INPUT:
% data       tracking info structure: cell mask, nucleus label matrix, and supporting information
% image_in    original cell image
% p           parameters structure from SetupTracking.m
%
% OUTPUT:
% output      all information (masks) needed for tracking
% diagnos     major masks/thresholds created as intermediates
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
diagnos = struct;
output.img_straight = image_in;
output.cells = data.nuclei;

% Save all information under diagnostic struct
diagnos = combinestructures(diagnos,output);