function [output, diagnos] =  nucleusID(nuc_orig,p,data)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [output, diagnos] =  nucleusID(nuc_orig,p,data,~) 
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% NUCLEUSID  Find nuclei from images of nuclear-localized fluorophore. Creates separated mask of identified nuclei.
% 
% nuc_orig        input fluorescent image
% p              parameters struture
% data           contains final cell mask from phaseID/ dicID (mask_cell)
%
% label_final    output mask showing cells 
% diag           structure with all masks and label matricies
%
%
% Subfunctions
% watershedalt.m, removemarked.m, bridgenuclei.m
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%- - - - - - - - - - - - - - - - - - - SETUP - - - - - - - - - - - - - - - - - - - - - - -
% Set cutoffs for nuclear shape
cutoff.Area = [floor(pi*(p.MinNucleusRadius)^2) ceil(pi*(p.MaxNucleusRadius)^2)];
cutoff.Compactness = p.Compactness;
cutoff.Solidity = p.Solidity;

if max(data.mask_cell)==0 % No cells found- break out of this
    output.label_nuc = zeros(size(data.mask_cell));
    diagnos = output;
    return;
end

% Pull out existing mask of cells
cell_mask = data.mask_cell;
% Add any strong nuclei (in case they weren't included in cell mask)
diagnos.thresh1 = quickthresh(nuc_orig,~cell_mask,'none');
tmp = nuc_orig>diagnos.thresh1;
if sum(tmp(:)) < sum(cell_mask(:))
    cell_mask = ~bwareaopen(~(cell_mask|tmp),p.NoiseSize,4);
end
% Construct smoothed images + watershed image
nucleus1 = double(medfilt2(uint16(nuc_orig),[p.MedianFilterSize, p.MedianFilterSize])); % Median-filtered

diagnos.nucleus_smooth1 = imfilter(nucleus1,gauss2D(p.MinNucleusRadius/4),'symmetric'); % Gaussian filtered
diagnos.watershed1 = watershedalt(diagnos.nucleus_smooth1, cell_mask, 4);
%- - - - - - - - - - - - - - - - - - - LABEL1: strong edges  - - - - - - - - - - - - - - - - - - - - - - -
% 1) Iterate down to p.NucleusEdgeThreshold to find strong-edge nuclei
horizontalEdge = imfilter(nucleus1,fspecial('sobel') /8,'symmetric');
verticalEdge = imfilter(nucleus1,fspecial('sobel')'/8,'symmetric');
diagnos.edge_mag = sqrt(horizontalEdge.^2 + verticalEdge.^2);
diagnos.edge_mag(nucleus1==max(nucleus1(:))) = max(diagnos.edge_mag(:)); % Correct for saturated nuclear centers
%%
tmp1 = diagnos.edge_mag(cell_mask);
edge_cutoffs = exp(linspace(log(p.NucleusEdgeThreshold), log(prctile(tmp1(:),97)),15));
cc_list = {};
for i = 1:length(edge_cutoffs)
    % a) Threshold, drop already-found objects
    mask0  = cell_mask & diagnos.edge_mag>=edge_cutoffs(end-i+1);
    tmp_drop = cell2mat(cc_list');
    if ~isempty(tmp_drop)
        mask0(tmp_drop) = 0;
    end
    % b) Skeletonize/ fill holes
    mask0 = bwmorph(mask0,'skel',2);
    mask0 = bwareaopen(mask0,p.NoiseSize,8);
    fill_size = cutoff.Area(2);
    if i >= (length(edge_cutoffs)-1)
        fill_size = round(6*cutoff.Area(1));
    end
    mask0 = ~bwareaopen(~mask0,fill_size,4); 
    if ~isempty(tmp_drop)
        mask0(tmp_drop) = 0;
    end
    % c) Filter objects that aren't round/sufficently large (alternate btw strict/lenient criteria)
    mask0 = imopen(mask0,diskstrel(round(p.NuclearSmooth)));
    if (mod(i-1,3) == 1) || (i==length(edge_cutoffs))
        mask0 = bwareaopen(mask0,round(cutoff.Area(1)),4);
    else
        mask0 = bwareaopen(mask0,2*cutoff.Area(1),4);
    end
    % d) Add newly-found objects to list
    cc_new = bwconncomp(mask0,8);
    cc_list = cat(2,cc_list,cc_new.PixelIdxList);
end

cc_all.PixelIdxList = cc_list';
cc_all.ImageSize = size(diagnos.edge_mag);
cc_all.NumObjects = length(cc_list);
cc_all.Connectivity = 4;
diagnos.label1a = labelmatrix(cc_all); % Edge-based division lines

% Quality check #1: ensure no object is wholly surrounded by another object, with no intervening background
erode1 = imerode(diagnos.label1a,ones(3));
dilate1 = imdilate(diagnos.label1a,ones(3));
get_unique = @(pixlist) unique([erode1(pixlist(:));dilate1(pixlist(:))]) ;
uniquelist = cellfun(get_unique,cc_all.PixelIdxList,'UniformOutput',0);
filter_obj = @(uniquelist) ~ismember(0,uniquelist) & (length(uniquelist)==2);
surroundeds = find(cellfun(filter_obj, uniquelist));
for i = 1:length(surroundeds)
    tmp_list = uniquelist{surroundeds(i)};
    diagnos.label1a(diagnos.label1a==surroundeds(i)) = tmp_list(tmp_list~=surroundeds(i));
end

%% 2) Label1b: subdivide objects using concave points on perimeter (~ >225 degrees)
tmp_label  =diagnos.label1a; tmp_label(diagnos.label1a==0) = max(diagnos.label1a(:))+1;
diagnos.mask_split = diagnos.label1a>0;
diagnos.mask_split((diagnos.label1a>0) & (tmp_label-imerode(tmp_label,ones(3)))>0)=0;
[mask_cut, diagnos.cut_pts] = perimetersplit(diagnos.mask_split,p);
diagnos.mask_split = diagnos.mask_split &~mask_cut;
diagnos.mask_split = bwareaopen(diagnos.mask_split,cutoff.Area(1),4);
diagnos.label1b = bwlabel(diagnos.mask_split,4);

%% 3) Label1c: subdivide objects with additional borders from edge-transformed image
if length(unique(diagnos.label1b(:)))>1
    
    if p.MinNucleusRadius < 6
       [diagnos.label1c, cc_obj] = edgesplit(diagnos.label1b, diagnos.watershed1, nucleus1, p);
    else
        [diagnos.label1c, cc_obj] = edgesplit(diagnos.label1b, diagnos.edge_mag, diagnos.cut_pts, p, mask_cut);
    end
    %diagnos.edge_borders = (border_mask) + (label_subobj>0);

    % Bridge oversegmented nuclear subobjects (from edge-based divisions) together by shape
    diagnos.label1 = bridgenuclei(diagnos.label1c, cc_obj, cutoff,p.ShapeDef, p.debug);
    
else % No objects were found- skip all these steps.
    diagnos.label1= diagnos.label1b;
end
%% - - - - - - - - - - - - - - - - - - - Label2 - - - - - - - - - - - - - - - - - - - - - - -
% ("Weak" objects missed by standard methods)
if p.WeakObjectCutoff>0
    % Drop mask1 "marked" watershed areas from watershed of Gaussian-smoothed image 
    label_dropped = imdilate(diagnos.watershed1,ones(3));
    markers = diagnos.label1>0;
    label_dropped = removemarked(label_dropped,markers,'remove')>0;
    nucleus_smooth2 = imfilter(nucleus1,gauss2D(p.MinNucleusRadius/2),'replicate'); % Gaussian filtered
    diagnos.watershed_remainder = imdilate(watershedalt(nucleus_smooth2, cell_mask, 4),ones(3));
    diagnos.watershed_remainder(label_dropped==0) = 0;

    % Rank remaining pixels, and use highest-valued pixels to bridge adjacent watershed regions
    diagnos.weak_ranked = rankpixels(diagnos.watershed_remainder, nucleus1);
    tmp_mask = diagnos.weak_ranked==max(diagnos.weak_ranked(:));
    tmp_mask = bwareaopen(tmp_mask,3); % Remove any speckle noise
    high_valued = bwconncomp(tmp_mask);
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
    % Remove small objects from (merged) watershed
    areas1 = cell2mat(struct2cell(regionprops(diagnos.watershed_remainder,'Area')));
    diagnos.watershed_remainder(ismember(diagnos.watershed_remainder,find(areas1<cutoff.Area(1)))) = 0;
    diagnos.weak_ranked2 = rankpixels(diagnos.watershed_remainder, nucleus1); % Rerank in merged watershed
    diagnos.watershed_remainder((imdilate(diagnos.watershed_remainder,ones(3))-diagnos.watershed_remainder)>0) = 0;

    % Check that brightest part of "nucleus" is relatively concentric-shaped and contiguous
    tmp_mask = imdilate(diagnos.weak_ranked2==1,ones(3))|(diagnos.weak_ranked==0);
    tmp_mask2 = diagnos.watershed_remainder==0;
    tmp_mask2(~tmp_mask) = 0;
    diagnos.weak_ranked2(tmp_mask2) = 0;
    test_weak =  diagnos.weak_ranked2 - imerode(diagnos.weak_ranked2,ones(3));

    % Label2a: based on watershed remainder
    diagnos.label2a = diagnos.watershed_remainder;
    diagnos.weak_objects = zeros(size(test_weak)); % (diagnostic image)
    weak_obj = label2cc(diagnos.label2a);
    
    get_score = @(pix) (sum(test_weak(pix)==2) + 3.3*sum(test_weak(pix)==3) + 10*sum(test_weak(pix)==4))/sqrt(numel(pix));
    all_scores = cellfun(get_score,weak_obj.PixelIdxList)/p.MinNucleusRadius; % scale measurement in range w/ prior score
    for i = 1:weak_obj.NumObjects
        diagnos.weak_objects(weak_obj.PixelIdxList{i}) = all_scores(i);
    end
    diagnos.weak_objects(diagnos.weak_objects>5) = 5;
    
    % Omit non-concentric objects, then get brightest 12% of pixels in region and proceed).
    diagnos.label2a(diagnos.weak_objects>p.WeakObjectCutoff) = 0;
    diagnos.label2a(diagnos.weak_ranked2<=2) = 0; 
    diagnos.label2a = imclose(diagnos.label2a,diskstrel(2));
    
    % Filter out small objects (label2b)
    diagnos.label2b = diagnos.label2a;
    diagnos.label2b(~imopen(diagnos.label2a>0,diskstrel(p.NuclearSmooth))) = 0;
    diagnos.label2b(~bwareaopen(diagnos.label2b,cutoff.Area(1))) = 0;
    % (Fix bug where some edge pixels belong to another object)
    diagnos.label2b= imerode(imdilate(diagnos.label2b,ones(3)),ones(3));
    diagnos.label2b= imdilate(imerode(diagnos.label2b,ones(3)),ones(3));
    diagnos.label2b = labelmatrix(label2cc(diagnos.label2b));
    cutoff.Area(1) = cutoff.Area(1)*0.5;
    
    % Do watershed + recombine using morphological properties
    diagnos.label2 = bridgenuclei(diagnos.label2b,bwconncomp(diagnos.label2b>0,4),cutoff,p.ShapeDef, p.debug);
else
    diagnos.label2 = zeros(size(diagnos.label1));
end

% - - - - - - - - - - - - - - - LABEL_END - - - - - - - - - - - - - - - - - - - - -
% Combine label1 and label2.
diagnos.label2(diagnos.label1>0) = 0; % Double check and make sure there's no overlap
diagnos.label2(diagnos.label2>0) = diagnos.label2(diagnos.label2>0)+max(diagnos.label1(:));
output.label_nuc = diagnos.label1+diagnos.label2;

% Relabel contiguously, just in case - convert to double so there are no math errors down the line,
tmp_cc = label2cc(output.label_nuc,1);
output.label_nuc = double(labelmatrix(tmp_cc));


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
    input_objects = label2cc(input_objects);
end

ranked_image = zeros(size(source_image));
for i = 1:input_objects.NumObjects
    locs = input_objects.PixelIdxList{i};
    if length(locs) < 20
        ranked_image(locs) = 4;
    else
        [~,sort_order] = sort(source_image(locs),'descend'); 
        vals = cat(2,4*ones(1,floor(0.08*length(locs))),3*ones(1,floor(0.12*length(locs))),2*ones(1,floor(0.15*length(locs))));
        vals = cat(2,vals,ones(1,length(locs)-length(vals)));
        ranked_image(locs(sort_order)) = vals;
    end
end


