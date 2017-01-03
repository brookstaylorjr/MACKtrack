function [output, diagnos] =  fluorescenceID(image_cell, p, image_nuc)
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

% Perform flatfield correction on images (if specified)
if p.CellFF>0    
    diagnos.img_cell = flatfieldcorrect(image_cell,double(p.Flatfield{p.CellFF}));
else
    diagnos.img_cell = image_cell;
end
if p.NucleusFF>0    
    diagnos.img_nuc = flatfieldcorrect(image_nuc,double(p.Flatfield{p.NucleusFF}));
else
    diagnos.img_nuc = image_nuc;
end

% Smooth and log-compress cell fluorescence image
img = imfilter(log(diagnos.img_cell+1),gauss2D(2),'symmetric');
diagnos.log_cell = img;
% Find image mode; determine whether it is a background or foreground intensity value
v = img(:);
bins = linspace(prctile(v,0.5),prctile(v,99.5),512);
N = histcounts(v,bins);
thresh1 = bins(find(N==max(N),1,'first'));
bin_width = bins(2)-bins(1);
contour_func = @(w) abs(sum(sum( (img<(thresh1 + bin_width*w))&(img>(thresh1- bin_width*w))))/numel(img(:)) - 0.05);
w = fminbnd(contour_func,0,20);
diagnos.mode_contours = (img<(thresh1 + bin_width*w)) &(img>(thresh1- bin_width*w));
contours_close = imclose(diagnos.mode_contours,diskstrel(4));
contours_open = imopen(contours_close,ones(2));
lost1 = (-sum(sum(contours_open))+sum(sum(contours_close)))/sum(sum(contours_close));

% Perform 1-D thresholding (stricter threshold; best if image mode is background)
diagnos.mask_nuc = diagnos.img_nuc  > quickthresh(diagnos.img_nuc,false(size(diagnos.img_nuc)),'log');
diagnos.mask0 = (img>trithreshold(img,bins,thresh1)) | diagnos.mask_nuc;

% If mode was determined to be foreground: use object-based (lenient) threshold 
if lost1> 0.17 % [threshold value separating contoured objects vs background - empirically determined, but may require adjusting]
    threshes= linspace(prctile(img(:),0.01),thresh1, 96);
        numobj = zeros(size(threshes));
        for j = 1:length(threshes)
            tmp = bwconncomp(img>threshes(j),4);
            numobj(j) = tmp.NumObjects;
        end
        v1 = 3:length(threshes)-3;
        tot_err = zeros(size(v1));
        for j = 1:length(v1)
            x1 = threshes(1:v1(j)); y1 = numobj(1:v1(j));
            f1 = polyfit(x1,y1,1);
            err1 = sum(abs(polyval(f1,x1)-y1));
            x2 = threshes(v1(j):end); y2 = numobj(v1(j):end);
            f2 = polyfit(x2,y2,1);
            err2 = sum(abs(polyval(f2,x2)-y2));
            tot_err(j) = err1+err2;
        end        
        diagnos.mask1 = (img>threshes(v1(find(tot_err==min(tot_err),1,'last')))) | diagnos.mask_nuc;
else
    diagnos.mask1 = diagnos.mask0;
end

% Morphological cleanup: remove speckle noise; do small close operation, and fill small holes
diagnos.mask_clean = bwareaopen(diagnos.mask1,2);
diagnos.mask_clean = imclose(diagnos.mask_clean,diskstrel(2));
diagnos.mask_clean = bwareaopen(diagnos.mask_clean,p.NoiseSize);
output.mask_cell = ~bwareaopen(~diagnos.mask_clean,p.MinHoleSize);

% (Output formatting)
diagnos = combinestructures(diagnos,output);
output.mask0 = diagnos.mask0; % Conservative estimate of cell boundaries - used in memory checking

function thresh = trithreshold(img,bins,mode_thresh)
v = img(:);

all_thresh = zeros(3,1);
% 1) Otsu threshold
all_thresh(1) = quickthresh(img,false(size(img)),'none');

% 2) Tsai threshold (based on histogram shape)
all_thresh(2) = tsaithresh(v, false(size(v)),length(bins));

% 3) 2-gaussian mix
start_2p = struct;
sigma1 = std(v(v<mode_thresh))*5/3; % Get std from bottom half of left distribution
start_2p.mu = [mode_thresh; mode_thresh+3*sqrt(sigma1)];
start_2p.Sigma = cat(3,sigma1, (sigma1*3));
p = min([0.97 numel(v(v<mode_thresh))*2/numel(v)]);
start_2p.PComponents = [p 1-p];
% Set maxIter for GM modeling (should be small for speed, <50)
gmopt = statset('MaxIter',50);
% Model distribution - subsample image a little bit (100K points)
dist_2p = fitgmdist(v(1:4:end),2,'Start',start_2p,'CovType','diagonal','Options',gmopt);
tmp = dist_2p.cluster(bins');
tmp_switch = [1+find(abs(diff(tmp))>0)'];
if length(tmp_switch)>1
    len1 = [tmp_switch(1)-1 length(tmp)-tmp_switch(2)];
    tmp_switch = tmp_switch(find(len1==max(len1),1,'first'));
end
all_thresh(3) = bins(tmp_switch);
thresh = median(all_thresh);