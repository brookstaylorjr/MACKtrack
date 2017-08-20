function [corrected_img, mult] = flatfieldcorrect(orig_img, flatfield)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% corrected_img = flatfieldcorrect(orig_img, flatfield)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% FLATFIELDCORRECT determines optimal scaling of a flatfield image (to reduce large-scale
% variation across a fluorescent image) and subtracts it from that image to correct it.
%
% INPUTS:
% orig_img   image to be corrected (should be double)
% flatfield  flatfield image- should be converted to 2-D cubic function using backgroundcalculate
%
% OUTPUTS:
% corr_img   flatfield-subtracted image
% mult       scaling factor for flatfield function
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

n = 4; % # of blocks ([n x n]) that image will be broken down into

blk_sz = [ceil(size(orig_img,1)/n) , ceil(size(orig_img,2)/n)];
lo_find = @(block_struct) prctile(block_struct.data(:),1);
guess = range(reshape(blockproc(orig_img,blk_sz,lo_find),[1 n^2]))...
    /range(reshape(blockproc(flatfield,blk_sz,lo_find),[1 n^2]));
% Calculating optimal scaling and subtract
fun1 = @(x) std(reshape(blockproc(orig_img./flatfield*x,blk_sz,lo_find),[1 n^2]));
opt1 = optimset('TolX',1e-2);
mult = fminbnd(fun1,0, guess*5,opt1);
corrected_img = orig_img ./ flatfield*mult;

% Rescale the corrected image to be 
