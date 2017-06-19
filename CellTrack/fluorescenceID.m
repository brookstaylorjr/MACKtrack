function [output, diagnos] =  fluorescenceID(image_cell, p, image_nuc)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [output, diagnos] =  fluorescenceID(image_cell, p, image_nuc) 
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
    if min(diagnos.img_cell(:))<0
        diagnos.img_cell = diagnos.img_cell-min(diagnos.img_cell(:));
    end
else
    diagnos.img_cell = image_cell;
end

% Mask nuclear image if it is provided, and selected for use.
if (nargin>2) && ~isempty(image_nuc)
    if p.NucleusFF>0    
        diagnos.img_nuc = flatfieldcorrect(image_nuc,double(p.Flatfield{p.NucleusFF}));
        if min(diagnos.img_nuc(:))<0
            diagnos.img_nuc = diagnos.img_nuc-min(diagnos.img_nuc(:));
        end
        diagnos.mask_nuc = diagnos.img_nuc > trithreshold(diagnos.img_nuc);
    else
        diagnos.mask_nuc = false(size(image_nuc));
    end
else
        diagnos.mask_nuc = false(size(image_cell));
end

% Smooth and log-compress cell fluorescence image
img = imfilter(log(diagnos.img_cell+1),gauss2D(2),'symmetric');
diagnos.log_cell = img;

% Get img histogram & mode
v = img(:);
try
    bins = linspace(prctile(v,0.5),prctile(v,99.5),512);
    N = histcounts(v,bins);
    N(1) = 0; N(end) = 0;
    thresh1 = bins(find(N==max(N),1,'first'));
    if p.Confluence == 2
        % Find image mode; determine whether it is a background or foreground intensity value
        bin_width = bins(2)-bins(1);
        contour_func = @(w) abs(sum(sum( (img<(thresh1 + bin_width*w))&(img>(thresh1- bin_width*w))))/numel(img(:)) - 0.05);
        w = fminbnd(contour_func,0,20);
        diagnos.mode_contours = (img<(thresh1 + bin_width*w)) &(img>(thresh1- bin_width*w));
        contours_close = imclose(diagnos.mode_contours,diskstrel(4));
        contours_open = imopen(contours_close,ones(2));
        lost1 = (-sum(sum(contours_open))+sum(sum(contours_close)))/sum(sum(contours_close));
        confluent = lost1 > 0.17;
    else
        confluent = p.Confluence;
    end

    % Perform 1-D thresholding (stricter threshold; best if image mode is background)
    if p.Confluence ~= 1 
        diagnos.mask0 = (img>trithreshold(img,bins,thresh1));
    end

    % If mode was determined to be foreground: use object-based (lenient) threshold 
    if confluent % [threshold value separating contoured objects vs background - empirically determined, but may require adjusting]
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
        diagnos.mask1 = (img>threshes(v1(find(tot_err==min(tot_err),1,'last'))));
    else
        diagnos.mask1 = diagnos.mask0;
    end

    % Add nuclear image to cell mask, if specified
    diagnos.mask2 = diagnos.mask1|diagnos.mask_nuc;
    
    
    % Morphological cleanup: remove speckle noise; do small close operation, and fill small holes
    diagnos.mask_clean = bwareaopen(diagnos.mask2,2);
    diagnos.mask_clean = imclose(diagnos.mask_clean,diskstrel(2));
    diagnos.mask_clean = bwareaopen(diagnos.mask_clean,p.NoiseSize);
    output.mask_cell = ~bwareaopen(~diagnos.mask_clean,p.MinHoleSize);

    % (Output formatting)
    diagnos = combinestructures(diagnos,output);

    if p.Confluence ~= 1
        output.mask0 = diagnos.mask0; % Conservative estimate of cell boundaries - used in memory checking
    else
        output.mask0 = diagnos.mask1;
    end
catch me % Fail gracefully (set all masks to false)
    output.mask_cell = false(size(image_cell));
    output.mask0 = false(size(image_cell));
    diagnos = combinestructures(diagnos,output);
end