
function [image_jump, r_jumps, c_jumps, maxes] = calculatejump(old_img, new_img, n)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [image_jump] = calculatejump(old_img, new_img)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% CALCULATEJUMP uses the normxcorr2 funtion to perform image registration between subsequent
% frames, calculating "jumps" between them. Offset calculation is given as shift from old
%  frame to new one
%
% INPUTS:
% old_img         frame (n-1) - prior to jump
% new_img         frame (n) - post-jump
% n               number of row/columns to divide image into - if not specified, defaults to 3 (a 3x3 grid)
%
% OUTPUT:
% image_jump      best-fit translation (rows, columns) from old frame to new frame
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%% Downsample input images - extra resolution beyond ~ 512x512 is unnecessary
if numel(old_img)> 2e6
    sz_down = floor(sqrt(numel(old_img) / 1e6));
    old_img = imresize(old_img,1/sz_down);
    new_img = imresize(new_img,1/sz_down);
else
    sz_down = 1 ;
end

%% Correlate blocks of new image w/ blocks of old one.
if nargin<3
    n = 3;
end

if length(n)<2
    n = [n n];
end

blocksize = ceil(size(new_img)./[n(1) n(2)]);
corr_fcn = @(block_struct) normxcorr2(block_struct.data(:,:,2),block_struct.data(:,:,1));
corr_2D = blockproc(cat(3,new_img,old_img),blocksize,corr_fcn,'PadPartialBlocks',1);

% Find maxima + correponding location for each block
blocksize2 = ceil(size(corr_2D)./[n(1) n(2)]);
get_max = @(block_struct) max(block_struct.data(:));
maxes = blockproc(corr_2D,blocksize2,get_max);
max_loc = @(block_struct) find(block_struct.data==max(block_struct.data(:)));
locs = blockproc(corr_2D,blocksize2,max_loc);
[r, c] = ind2sub(blocksize2,locs);


% Subtract reference points to compute overall jump
r_jumps = (r - blocksize(1))*sz_down;
c_jumps = (c - blocksize(2))*sz_down;


all_jumps = [r_jumps(:) c_jumps(:)];
all_jumps(maxes(:)<0.35,:) = nan;  
image_jump = nanmedian(all_jumps);

image_jump(isnan(image_jump)) = 0; % If we couldn't get an accurate fix, just assume no jump at all.






