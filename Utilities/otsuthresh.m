function lvl = otsuthresh(inImg,dropPixels,logComp)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% OTSUTHRESH Perform Otsu thresholding on image using modified form of 'graythresh' function
%
%  [lvl] = otsuthresh(inImg,dropPixels,logComp)
%
% inImg        the input image to be thresholded
% dropPixels   any image pixels that are not included in the calculation
% logComp      'none' (no modification) or 'log' (image is log-compressed)
%
% lvl   calculated threshold of image
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

im= inImg(:);

if nargin<3
    logComp = 'none';
    if nargin<2
        dropPixels = false(size(im));
    end
end
   
im(dropPixels(:)) = [];
        
if strcmp(logComp,'log')
    im(im<=0) = min(im>0);
    im = log(im);
    minval = min (im);
    maxval = max (im);
    im = (im - minval) / (maxval - minval);
    lvl = exp(minval + (maxval - minval) * graythresh(im));
else
    minval = min (im);
    maxval = max (im);
    im = (im - minval) / (maxval - minval);
    lvl = (minval + (maxval - minval) * graythresh(im));
end

