function [output, diagnos] =  fluorescenceID(image_cell,p, image_nuc)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [output, diagnos] =  fluorescenceID(image0,p,data,~) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% FLUORESCENCEID  Identifies candidate regions for cells (foreground) in a fluorescent image
%
% INPUTS:
% image0        input fluorescent (cell) image
% p             parameters struture
% image1        input fluorescent (nucleus) image

%
% OUTPUTS:
% output        all masks needed for tracking purposes
% diagnos       structure with all masks and label matricies (for diagnosis in testing mode)
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%- - - - - - - - - - - - - - - - - - - MAKE "CELL" MASK - - - - - - - - - - - - - - - - - - - - - - -
%%
% Normalize images via flatfield correction -> estimate range for scaling
img = {image_cell, image_nuc};
corrected = cell(2,1);
for i = 1:length(img)
    orig_img = double(img{i}); bg_img = double(p.Flatfield{i});
    blk_sz = [ceil(size(orig_img,1)/6) , ceil(size(orig_img,2)/6)];
    lo_find = @(block_struct) prctile(block_struct.data(:),2);
    guess = range(reshape(blockproc(orig_img,blk_sz,lo_find),[1 36]))/range(reshape(blockproc(bg_img,blk_sz,lo_find),[1 36]));
    % Calculate optimal scaling and subtract
    fun1 = @(x) std(reshape(blockproc(orig_img-bg_img*x,blk_sz,lo_find),[1 36]));
    opt1 = optimset('TolX',1e-4);
    mult = fminbnd(fun1,0, guess*5,opt1);
    corrected{i} = orig_img - bg_img*mult;
end

diagnos.cell_corr = corrected{1};
diagnos.nuc_corr = corrected{2};



% Get (very) conservative threshold and use to define cell mask
[~,bg_dist] = modebalance(diagnos.cell_corr, 3, 16, 'measure');
thresh1 = min([ quickthresh(diagnos.cell_corr+10,false(size(diagnos.cell_corr)),'log'), bg_dist(1)]);
diagnos.mask1 = diagnos.cell_corr>thresh1;

% Add in nuclear mask (copies primaryID)
[~,bg_dist] = modebalance(diagnos.nuc_corr, 2, 16, 'measure');
thresh1 = min([ quickthresh(diagnos.nuc_corr,false(size(diagnos.nuc_corr)),'none'), bg_dist(1)+4*bg_dist(2)]);
diagnos.mask2 = diagnos.nuc_corr>thresh1;



% Morphological cleanup: speckle remove, closing, larger subobject removal, then cleanup.
diagnos.mask_clean = bwareaopen(diagnos.mask1|diagnos.mask1,2);
diagnos.mask_clean = imclose(diagnos.mask_clean,diskstrel(3));
diagnos.mask_clean = bwareaopen(diagnos.mask_clean,p.NoiseSize);

diagnos.mask_fill = ~bwareaopen(~diagnos.mask_clean,p.MinHoleSize);

output.mask0 = diagnos.mask_clean;
output.mask_cell = diagnos.mask_fillfi;


