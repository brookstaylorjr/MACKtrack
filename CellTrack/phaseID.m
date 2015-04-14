function [output, diagnos] =  phaseID(phaseOrig,p,X)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% PHASEID:   Edge and intensity method of determining background vs cell.
%
% phaseOrig    original phase-contrast image
% p            parameters struture
% X            background correction matix multiplier 
%
% output       masks used in downstream tracking
% diagnos      structure with all intermediate masks and thresholds
%
% Structure: We operate on two types of data. The original image is background corrected, then thresholded, showing
% "halo" artifacts from the image. The image is also edge-transformed using Sobel horizontal/vertical filters, smoothed
% slightly, and thresholded. This mask is used as a rough identifier of foreground (cells): closing and hole filling
% refine the mask, and it is "walked in" from the outside, using both edge datarmation and data.halos, yielding a final mask.
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Subfunctions
% otsuthresh.m, tsaithresh.m, noisethresh.m, fillholes.m
%

% Notes
% 12/14/2010 - Created from edgeID-Hasty.m. Remaining non-parametized values are:
%       - imgSubset dilation size (set at 80)
%       - speckle noise size (4)
%       - small object size (512)
% 01/19/2011- Converted script into function, made modifications to hole filling and edge determination 
% 02/16/2011- Changed initial edge determination method, using marking of Canny edges
% 03/23/2011- Modified gap-filling and hole-filling method
% 01/09/2012- Changed output variable names

%__________________________ 1. Load and background-correct______________________________________
% Calculate background function
warning off MATLAB:nearlySingularMatrix
pStar = (X'*X)\(X')*double(phaseOrig(:));
% Apply correction
diagnos.phase1 = reshape((double(phaseOrig(:) - X*pStar)),size(phaseOrig));
diagnos.phase1 = diagnos.phase1-min(diagnos.phase1(:)); % Set minimum to zero

%_______________________ 2. Edge transformation on image_________________________________
horizontalEdge = imfilter(phaseOrig,fspecial('sobel') /8,'replicate');
verticalEdge = imfilter(phaseOrig,fspecial('sobel')'/8,'replicate');
diagnos.combined_edge = sqrt(horizontalEdge.^2 + verticalEdge.^2);

%________________________ 3. Control for sparsely populated images ___________________________
diagnos.subsetThreshold = otsuthresh(diagnos.combined_edge,false(size(diagnos.combined_edge)),'none');
diagnos.image_subset = imdilate(diagnos.combined_edge>(diagnos.subsetThreshold), ones(80));


% _______________________ 4. Intensity/edge thresholding _______________________________________
diagnos.haloThreshold = tsaithresh(diagnos.phase1,~diagnos.image_subset);
output.halos = diagnos.phase1>(diagnos.haloThreshold);
[diagnos.edgeThreshold diagnos.noiseCount diagnos.val]  = noisethresh(diagnos.combined_edge, ~diagnos.image_subset, p.PhaseSearchRange ,p.NoiseSize);
diagnos.mask_orig = diagnos.combined_edge>diagnos.edgeThreshold;
diagnos.mask_orig = bwareaopen(diagnos.mask_orig,3); % Get rid of outlying noise

% _______________________5. Masking, hole and gap-filling ___________________________________
% Use strong edges to mark/keep Canny edges (7-pixel filter)
thresh1 = [0.0125    0.0312];
BWcannyFilt7 = edge(diagnos.phase1,'canny',thresh1,3);
BWcannyFilt7 = imdilate(BWcannyFilt7,ones(3));
BWcannyFilt7 = bwmorph(BWcannyFilt7,'skel','Inf');
edgeMarkers = bwareaopen(BWcannyFilt7 & diagnos.mask_orig,2);
diagnos.mask_canny7 = diagnos.mask_orig | (BWcannyFilt7 &~ (removemarked(bwlabel(BWcannyFilt7),edgeMarkers)>0));

% Fill holes in connected mask, being lenient on small bright holes
diagnos.mask_filled1 = fillholes(diagnos.mask_canny7,output.halos,p.HaloCutoff/100,p.MinHoleSize,p.MaxHoleSize);
diagnos.mask_filled1 = bwareaopen(diagnos.mask_filled1,p.NoiseSize);

% Final "connection" step: dilate out remaining edges, thin result (cycle through structuring elements)
% SEround = diskstrel(floor(p.MinCellWidth/2));
% BWremain = data.mask_filled1 &~ imopen(data.mask_filled1, SEround);
% BWremain = bwareaopen(BWremain,p.NoiseSize/2);
BWremain = bwmorph(diagnos.mask_filled1,'skel','Inf');
for i = 3:2:floor(sqrt(p.NoiseSize))
    BWremain = imdilate(BWremain,ones(i));
    BWremain = bwmorph(BWremain,'skel','Inf');
    % Reduce back result of dilation
    for j = 1:floor(i/2)
        BWremain(bwmorph(BWremain,'endpoints')) = 0;
    end
end

% Fill holes in this connected mask
diagnos.mask_filled2 = fillholes(BWremain|diagnos.mask_filled1,output.halos,p.HaloCutoff/100,0,p.MaxHoleSize);
diagnos.mask_filled2 = imclose(diagnos.mask_filled2,ones(3));


% _________________________7. Walk edges in _______________________________________________
% Define edges for walking in mask
lowEdge = prctile(abs([horizontalEdge(diagnos.image_subset);verticalEdge(diagnos.image_subset)]),50); % Use 50th percentile of combined vert/hor edge as threshold

diagnos.horizontal_test = -1*(bwareaopen(horizontalEdge<-lowEdge,4)) + (bwareaopen(horizontalEdge>lowEdge,4));
diagnos.vertical_test = -1*(bwareaopen(verticalEdge<-lowEdge,4)) + (bwareaopen(verticalEdge>lowEdge,4));
[BWphase, diagnos.walked_in] = walkin(diagnos.mask_filled2,diagnos.vertical_test,diagnos.horizontal_test,output.halos, p.LeakThrough, p.WalkIn1,p.WalkIn2);

% Remove noise for final mask
output.mask_cell = bwareaopen(BWphase,p.NoiseSize,4);

% Save all information under diagnostic struct
diagnos = combinestructures(diagnos,output);

% ========================================================================================



function [outMask,diagnosticImg] = walkin(inMask,vertEdges,horEdges,intensityMask, leakThrough, iterations1,iterations2)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% WALKIN Take an input mask and iteratively shrink the boundaries by logic applied to top/bottom/left/right edges.  
%
% inMask         original mask
% vertEdges      smoothed vertical Sobel transformation
% horEdges       smoothed horiz. Sobel transformation
% intensityMask  intensity-thresholded image 
% leakThrough    parameter sets minimum eroded size
% iterations1    iterations from bottom/right
% iterations2    from top/left, sped up version of bottom/right
%
% outMask        new, shrunk version of original mask
% diagnosticImg  diagnostic image showing which pixels were removed in each iteration
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

diagnosticImg = zeros(size(inMask));
bottomEdges = false(size(inMask));
topEdges = false(size(inMask));
leftEdges = false(size(inMask));
rightEdges = false(size(inMask));
outMask = inMask;

for i = 1:iterations1+iterations2
    b = diff(outMask,1,2);
    a = diff(outMask,1,1);

    bottomEdges(1:end-1,1:end) = a<0;
    rightEdges(1:end,1:end-1) = b<0;
     topEdges(2:end,1:end) = a>0;
     leftEdges(1:end,2:end) = b>0;
    % General rule: we assume cell edges are shown by edge direction
    %                                 (+)
    %                      (+)         o        (-)
    %                                 (-)
    % On later iterations, we also pull out bright pixels (halos)
    erodeMat1= (topEdges &  (horEdges<=0) )...
                 | (leftEdges & (vertEdges<=0));
    erodeMat2 = (bottomEdges & (horEdges>=0))...
                | (rightEdges &  (vertEdges>=0));

    if i >iterations1    
        erodeMat1= (topEdges &  (intensityMask |(horEdges<=0) ))...
                 | (leftEdges & (intensityMask | (vertEdges<=0))) ;
        erodeMat2 = (bottomEdges & (intensityMask |(horEdges>=0)))...
                | (rightEdges &  (intensityMask | (vertEdges>=0))) ;
    end

    erodeMask2 = bwareaopen(erodeMat2,leakThrough); 
    erodeMask1 = bwareaopen(erodeMat1,leakThrough); 
   outMask = ~(erodeMask1|erodeMask2) & outMask;
    
    diagnosticImg = diagnosticImg+outMask;
end
% ========================================================================================


function BW1 = fillholes(BW0,halos,cutoff,minTotal,maxTotal)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% FILLHOLES:  Fills holes in mask as long as they do not contain too many "halo" pixels
%
% BW0       unfilled mask
% halos     thresholded intensity image
% cutoff    maximum proportion of "bright" pixels in a hole
% minTotal  fill all holes smaller than this value
% maxTotal  fill all holes larger than this value
%
% BW1       filled mask
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


BW1 = ~bwareaopen(~BW0,maxTotal,4);
allHoles = bwconncomp(BW1 & ~BW0,4);

for i = 1:allHoles.NumObjects
    if ((sum(halos(allHoles.PixelIdxList{i}))/length(allHoles.PixelIdxList{i})) > cutoff)&&(length(allHoles.PixelIdxList{i})>minTotal)
        BW1(allHoles.PixelIdxList{i}) = 0;
    end
end