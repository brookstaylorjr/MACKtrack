function [out_array] = smoothrows(in_array, smooth_sz)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [out_array] = smoothrows(in_array, smooth_sz)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SMOOTHROWS provides fast smoothing across rows of an array.
%
% in_array    array to be smoothed
% smooth_sz   size of smoothing window
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if nargin<2
    smooth_sz = 3;
end

tmp = [nan(size(in_array,1),ceil(smooth_sz/2)),in_array,nan(size(in_array,1),ceil(smooth_sz/2))]';
c = smooth(tmp(:),smooth_sz,'moving');
out_array = reshape(c,size(tmp,1),size(tmp,2))';
out_array = out_array(:,ceil(smooth_sz/2)+1:end-ceil(smooth_sz/2));