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


% Break new image up into subregions - try to cross-correlate with old image.

rows = round(linspace(1,size(new_img,1),9));
cols = round(linspace(1,size(new_img,2),9));
z = [];
all_diffs = nan(8,2);
for i = 1:8
    img_corner = [rows(i) cols(i)];
    template = new_img(rows(i):rows(i+1),cols(i):cols(i+1));
    tmp = normxcorr2(template,old_img);
    if max(tmp(:))>0.6
        [tmp_r, tmp_c] = find(tmp==max(tmp(:)));
        all_diffs(i,:) = img_corner-[tmp_r-size(template,1), tmp_c-size(template,2)];
    end
end
image_jump = round(nanmedian(all_diffs));
image_jump(isnan(image_jump)) = 0;