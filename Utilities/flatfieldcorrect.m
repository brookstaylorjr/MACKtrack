function [corrected_img, mult] = flatfieldcorrect(orig_img, flatfield,mode)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% corrected_img = flatfieldcorrect(orig_img, flatfield,mode)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% FLATFIELDCORRECT determines optimal scaling of a flatfield image (to reduce large-scale
% variation across a fluorescent image) and subtracts it from that image to correct it.
%
% INPUTS:
% orig_img   image to be corrected (should be double)
% flatfield  flatfield image- should be converted to 2-D cubic function using backgroundcalculate
% mode       'divide' or 'subtract' method of correction - 'divide' is default
%
% OUTPUTS:
% corr_img   flatfield-subtracted image
% mult       scaling factor for flatfield function
%
% Note on correction method: divide is more accurate overall, especially for cell measurements,
% but subtraction tends to work better in proccessing background pixels accurately - it may be
% preferred for object masking.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin<3
    mode = 'divide';
end
assert(sum(strcmp({'divide', 'subtract'},mode))>0,[ 'Error: correction mode (3rd argument) must be ', ...
    'either ''divide'' or ''subtract'''])

if strcmp(mode,'divide')
    corrected_img = orig_img ./ flatfield;
    corrected_img = corrected_img*mean(flatfield(:)); % Rescale image so it's approximately @ original intensity level

elseif strcmp(mode,'subtract')
    % In "subtract" mode, need to rescale flatfield 1st so that it's on same scale as acquired image.
    n = 6; % # of blocks ([n x n]) that image will be broken down into
    blk_sz = [ceil(size(orig_img,1)/n) , ceil(size(orig_img,2)/n)];
    lo_find = @(block_struct) prctile(block_struct.data(:),1);
    guess = range(reshape(blockproc(orig_img,blk_sz,lo_find),[1 n^2]))...
        /range(reshape(blockproc(flatfield,blk_sz,lo_find),[1 n^2]));
    fun1 = @(x) std(reshape(blockproc(orig_img - flatfield*x,blk_sz,lo_find),[1 n^2]));
    opt1 = optimset('TolX',1e-2);
    mult = fminbnd(fun1,0, guess*5,opt1);
    corrected_img = orig_img - flatfield*mult;
end

