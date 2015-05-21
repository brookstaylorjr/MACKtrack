function [output, diagnos] =  nucleusID(nucOrig,p,data,~)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% NUCLEUSID  Find nuclei from images of nuclear-localized fluorophore. Creates separated mask of identified nuclei.
% 
% nucOrig        input fluorescent image
% p              parameters struture
% data           contains final cell mask from phaseID/ dicID (mask_cell)
%
% label_final    output mask showing cells 
% diag           structure with all masks and label matricies
%
%
% Subfunctions
% watershedalt.m, removemarked.m
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%- - - - - - - - - - - - - - - - - - - SETUP - - - - - - - - - - - - - - - - - - - - - - -
% Set cutoffs for nuclear shape
cutoff.Area = [floor(pi*(p.MinNucleusRadius-1)^2) ceil(pi*(p.MaxNucleusRadius)^2)];
cutoff.Eccentricity = p.Eccentricity;
cutoff.Compactness = p.Compactness;

% Pull out existing mask of cells
cell_mask = data.mask_cell;
% Add any strong nuclei (in case they weren't included in cell mask)
diagnos.thresh1 = otsuthresh(nucOrig,~cell_mask,'none');
tmp = nucOrig>diagnos.thresh1;
if sum(tmp(:)) < sum(cell_mask(:))
    cell_mask = ~bwareaopen(~(cell_mask|tmp),p.NoiseSize,4);
end
% Construct smoothed images + watershed image
nucleus1 = medfilt2(nucOrig,[p.MedianFilterSize, p.MedianFilterSize]); % Median-filtered
diagnos.nucleus_smooth1 = imfilter(nucleus1,gauss2D(floor(p.MinNucleusRadius/4)),'replicate'); % Gaussian filtered
diagnos.watershed1 = watershedalt(diagnos.nucleus_smooth1, cell_mask, 4);

%- - - - - - - - - - - - - - - - - - - Label1 - - - - - - - - - - - - - - - - - - - - - - -
% Label1: strong-edge nuclei. Use watershed1 and p.NucleusEdgeThreshold
horizontalEdge = imfilter(nucleus1,fspecial('sobel') /8,'symmetric');
verticalEdge = imfilter(nucleus1,fspecial('sobel')'/8,'symmetric');
diagnos.edge_mag = sqrt(horizontalEdge.^2 + verticalEdge.^2);
diagnos.edge_mag(nucleus1==max(nucleus1(:))) = max(diagnos.edge_mag(:)); % Correct for saturated nuclear centers

% Take subset of edge vals > p.NucleusEdgeThresh- step down incrementally
edge_cutoffs = prctile(diagnos.edge_mag(diagnos.edge_mag>p.NucleusEdgeThreshold),0:20:80);
unique_all = 0;
diagnos.label1a = zeros(size(diagnos.watershed1));
for i = 1:length(edge_cutoffs)
    % Threshold, fill holes, check for minsize
    mask_tmp = ~bwareaopen(diagnos.edge_mag<edge_cutoffs((end-i)+1),cutoff.Area(2),4);
    mask_tmp = imopen(mask_tmp&cell_mask,diskstrel(p.MinNucleusRadius));
    unique_tmp = unique(diagnos.watershed1(mask_tmp));
    % If we found object before, drop it.
    unique_tmp(ismember(unique_tmp,unique_all)) = [];
    % Record new object positions, then add new unique objects to list.
    for j = 1:length(unique_tmp)
       diagnos.label1a((diagnos.watershed1==unique_tmp(j)) & mask_tmp) = unique_tmp(j);
    end
    unique_all = cat(1,unique_all, unique_tmp(:));
end

diagnos.label1a = labelmatrix(label2cc(diagnos.label1a));
diagnos.label1 = bridgenuclei(diagnos.label1a,cutoff,p.debug);

% ID nuclei w/ strong edges tends to be over-generous. Erode things somewhat, then remove super-small objects again
borders = (imdilate(diagnos.label1,ones(3))-diagnos.label1)>0;
strelsize = floor(p.MinNucleusRadius/3);
diagnos.label1(imdilate(borders,diskstrel(strelsize))) = 0;
diagnos.label1(~bwareaopen(diagnos.label1>0,cutoff.Area(1),4)) = 0;

%- - - - - - - - - - - - - - - - - - - Label2 - - - - - - - - - - - - - - - - - - - - - - -
% Drop mask1 "marked" watershed areas from watershed of Gaussian-smoothed image 
label_dropped = imdilate(diagnos.watershed1,ones(3));
markers = diagnos.label1>0;
label_dropped = removemarked(label_dropped,markers,'remove')>0;
nucleus_smooth2 = imfilter(nucleus1,gauss2D(p.MinNucleusRadius/2),'replicate'); % Gaussian filtered
diagnos.watershed_remainder = imdilate(watershedalt(nucleus_smooth2, cell_mask, 4),ones(3));
diagnos.watershed_remainder(label_dropped==0) = 0;

% Rank remaining pixels, and use highest-valued pixels to bridge adjacent watershed regions
diagnos.weak_ranked = rankpixels(diagnos.watershed_remainder, nucleus1);
high_valued = bwconncomp(diagnos.weak_ranked==max(diagnos.weak_ranked(:)));
% Merge watershed areas based on connected "high" areas
for i = 1:high_valued.NumObjects
    obj = unique(diagnos.watershed_remainder(high_valued.PixelIdxList{i}));
    obj(obj==0) = [];
    if length(obj)>1
        for j = 2:length(obj)
            diagnos.watershed_remainder(diagnos.watershed_remainder==obj(j)) = obj(1);
        end
    end
end
diagnos.watershed_remainder((imdilate(diagnos.watershed_remainder,ones(3))-diagnos.watershed_remainder)>0) = 0;
diagnos.weak_ranked2 = rankpixels(diagnos.watershed_remainder, nucleus1); % Rerank in merged watershed


% Check that brightest part of "nucleus" is relatively concentric-shaped and contiguous
test_weak =  diagnos.weak_ranked2 - imerode(diagnos.weak_ranked2,ones(3));
bright_edge = bwareaopen(test_weak==4,8);
test_weak(bright_edge) = 100; % Penalize cells with strong intensity values near edge


% Label2a: based on watershed remainder
diagnos.label2a = diagnos.watershed_remainder;
diagnos.weak_objects = zeros(size(test_weak)); % (diagnostic image)
weak_obj = label2cc(diagnos.label2a);
for i = 1:weak_obj.NumObjects
    testval = mean(test_weak(weak_obj.PixelIdxList{i}));
    diagnos.weak_objects(weak_obj.PixelIdxList{i}) = min([testval,3]);
    if (testval > p.WeakObjectCutoff) || (testval==0)
        diagnos.label2a(weak_obj.PixelIdxList{i}) = 0;
    end      
end
% Clean up label2
diagnos.label2a(diagnos.weak_ranked2<=2) = 0; % Only look at brightest 25% of area
diagnos.label2a = imclose(diagnos.label2a,diskstrel(2));
diagnos.label2a(~imopen(diagnos.label2a>0,diskstrel(floor(p.MinNucleusRadius*2/3)))) = 0;

% Fix bug where some edge pixels belong to another object
diagnos.label2a= imerode(imdilate(diagnos.label2a,ones(3)),ones(3));
diagnos.label2a= imdilate(imerode(diagnos.label2a,ones(3)),ones(3));
diagnos.label2a = labelmatrix(label2cc(diagnos.label2a));
cutoff.Area(1) = cutoff.Area(1)*0.5;
diagnos.label2 = bridgenuclei(diagnos.label2a,cutoff,p.debug);


% - - - - - - - - - - - - - - - LABEL_END - - - - - - - - - - - - - - - - - - - - -
% Combine label1 and label2.
diagnos.label2(diagnos.label1>0) = 0; % Double check and make sure there's no overlap
diagnos.label2(diagnos.label2>0) = diagnos.label2(diagnos.label2>0)+max(diagnos.label1(:));
output.label_nuc = diagnos.label1+diagnos.label2;

% Save all information under diagnostic struct
diagnos = combinestructures(diagnos,output);


%==================================================================================================



function ranked_image = rankpixels(input_objects, source_image)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% RANKPIXELS outputs a ranked image- input_objects (a bwconncomp structure) are sorted 
% into divisions- the highest 10% pixels are assigned "4", the next 15% are assigned "3",
% then "2", etc.
%
% input_objects       bwconncomp of objects
% source_image        matrix of values to rank
%
% ranked_image        equivalent to size of source_image, can take value of 0-3
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if ~isstruct(input_objects)
    labelmat = input_objects;
    all_obj = unique(labelmat);
    input_objects = struct;
    input_objects.NumObjects = length(all_obj)-1;
    input_objects.PixelIdxList = cell(1,length(all_obj)-1);
    for i = 1:(length(all_obj)-1)
        input_objects.PixelIdxList{i} = find(labelmat==all_obj(i+1));
    end
end

ranked_image = zeros(size(source_image));
for i = 1:input_objects.NumObjects
    locs = input_objects.PixelIdxList{i};
    if length(locs) < 20
        ranked_image(locs) = 4;
    else
        [~,sort_order] = sort(source_image(locs),'descend'); 
        vals = cat(2,4*ones(1,floor(0.1*length(locs))),3*ones(1,floor(0.15*length(locs))),2*ones(1,floor(0.2*length(locs))));
        vals = cat(2,vals,ones(1,length(locs)-length(vals)));
        ranked_image(locs(sort_order)) = vals;
    end
end


