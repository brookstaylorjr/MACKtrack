
function [image_jump] = calculatejump(old_img, new_img)
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
%
% OUTPUT:
% image_jump      best-fit translation (rows, columns) from old frame to new frame
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%% Downsample input images - extra resolution beyond ~ 512x512 is unnecessary
if numel(old_img)> 1e6
    sz_down = floor(sqrt(numel(old_img) / 2.5e5));
    old_img = imresize(old_img,1/sz_down);
    new_img = imresize(new_img,1/sz_down);
else
    sz_down = 1 ;
end
%%
% Correlate blocks of new image w/ old one
n = 3;
blocksize = ceil(size(new_img)/n);
corr_fcn = @(block_struct) normxcorr2(block_struct.data,old_img);
corr_2D = blockproc(new_img,blocksize,corr_fcn);

% Find maxima + correponding maxima for each block
blocksize2 = ceil(size(corr_2D)/n);
get_max = @(block_struct) max(block_struct.data(:));
maxes = blockproc(corr_2D,blocksize2,get_max);
max_loc = @(block_struct) find(block_struct.data==max(block_struct.data(:)));
locs = blockproc(corr_2D,blocksize2,max_loc);


% Get reference points for each block used in correlation
get_corners = @(block_struct) block_struct.location;
corners = blockproc(old_img,blocksize,get_corners);
[r, c] = ind2sub(blocksize2,locs);
r_corner = corners(:,1:2:end);
c_corner = corners(:,2:2:end);
    
image_jump = -[r(:) - r_corner(:)-blocksize(1), c(:) - c_corner(:)-blocksize(2)];
image_jump(maxes(:)<0.35) = nan;
image_jump = nanmedian(image_jump)*sz_down;
image_jump(isnan(image_jump)) = 0; % If we couldn't get an accurate fix, just assume no jump at all.






