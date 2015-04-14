function kk = lpfilter(sz, subsmp, n)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% LPFILTER is a simple lowpass filter kernel (Butterworth).
%
% sz      size of the filter
% subsmp  downsampling factor to be used later
% n       degree of the butterworth filter
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

sz = 2*floor(sz/2)+1; % make sure the size of the filter is odd
cut_frequency = 0.5 / subsmp;
range = (-(sz-1)/2:(sz-1)/2)/(sz-1);
[ii,jj] = ndgrid(range,range);
rr = sqrt(ii.^2+jj.^2);
kk = ifftshift(1./(1+(rr./cut_frequency).^(2*n)));
kk = fftshift(real(ifft2(kk)));
kk = kk./sum(kk(:));