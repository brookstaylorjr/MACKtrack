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

% Pull out existing mask of cells
cell_mask = data.mask_cell;
% Add any strong nuclei (in case they weren't included in cell mask)
diagnos.thresh1 = quickthresh(nuc_orig,~cell_mask,'none');
tmp = nuc_orig>diagnos.thresh1;
if sum(tmp(:)) < sum(cell_mask(:))
    cell_mask = ~bwareaopen(~(cell_mask|tmp),p.NoiseSize,4);
end
% Construct smoothed images + watershed image
nucleus1 = medfilt2(nuc_orig,[p.MedianFilterSize, p.MedianFilterSize]); % Median-filtered

diagnos.nucleus_smooth1 = imfilter(nucleus1,gauss2D(p.MinNucleusRadius/4),'replicate'); % Gaussian filtered
diagnos.watershed1 = watershedalt(diagnos.nucleus_smooth1, cell_mask, 4);
%- - - - - - - - - - - - - - - - - - - LABEL1: strong edges  - - - - - - - - - - - - - - - - - - - - - - -
% 1) Iterate down to p.NucleusEdgeThreshold to find strong-edge nuclei
horizontalEdge = imfilter(nucleus1,fspecial('sobel') /8,'symmetric');
verticalEdge = imfilter(nucleus1,fspecial('sobel')'/8,'symmetric');
diagnos.edge_mag = sqrt(horizontalEdge.^2 + verticalEdge.^2);
diagnos.edge_mag(nucleus1==max(nucleus1(:))) = max(diagnos.edge_mag(:)); % Correct for saturated nuclear centers
tmp1 = diagnos.edge_mag(cell_mask);
edge_cutoffs = linspace(p.NucleusEdgeThreshold, prctile(tmp1(:),95),21);
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
    mask0 = ~bwareaopen(~mask0,cutoff.Area(2)*4,4);    
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
    % Isolate high edge (>median) pixels per object
    tmp_cc = label2cc(diagnos.label1b,0);
    get_high = @(pix) pix(diagnos.edge_mag(pix) > prctile(diagnos.edge_mag(pix),50));
    mask_vals = cell2mat(cellfun(get_high,tmp_cc.PixelIdxList,'UniformOutput',0));
    mask_hi = false(size(diagnos.label1b));
    mask_hi(mask_vals) = 1;
    % Morphological cleanup and connection - dilate out, then shrink back result
    mask_hi = bwareaopen(mask_hi,round(p.MinNucleusRadius/2));
    mask_hi = ~bwareaopen(~mask_hi,round(p.MinNucleusRadius+1),4);
    BWconnect = bwmorph(mask_hi,'skel','Inf');
    for i = 3:2:round(sqrt(p.MinNucleusRadius))
        BWconnect = imdilate(BWconnect,ones(i));
        BWconnect = bwmorph(BWconnect,'skel','Inf');
        % Reduce back result of dilation
        for j = 1:floor(i/2)
            BWendpoints = bwmorph(bwmorph(BWconnect,'endpoints'),'shrink',Inf);
            BWconnect(BWendpoints) = 0;
        end
    end
    mask_hi = BWconnect;
    mask_obj = diagnos.label1b>0;
    off_mask = mask_hi|BWconnect;
    % For "solid" objects (i.e. entirely high-edge), get borders of these directly and substitute
    solid_mask = imopen(mask_hi,diskstrel(p.NuclearSmooth));
    solid_borders = solid_mask&~imerode(solid_mask,ones(3));
    off_mask(solid_mask) = solid_borders(solid_mask);
    % Break up existing mask w/ newly-found borders -> use propagate subfcn to fill existing mask.
    mask_obj(off_mask) = 0;
    mask_obj = imopen(mask_obj,ones(2));
    mask_obj = bwareaopen(mask_obj,round(cutoff.Area(1)/4),4);
    mask_all = diagnos.label1b>0;
    diagnos.label1c = IdentifySecPropagateSubfunction(double(bwlabel(mask_obj,4)),double(mask_all),mask_all,0.02);
    
    % Combine inflection point & strong edge data - see if a strong edge unabigiously connects two moderately-inflected pts.
    % (Filter out objects that are too small to be split further)
    filter_areas = @(pix) length(pix)<(2*cutoff.Area(1));
    nosmall = label2cc(diagnos.label1b);
    nosmall.PixelIdxList(cellfun(filter_areas,nosmall.PixelIdxList)) = [];
    nosmall.NumObjects = length(nosmall.PixelIdxList);
    nosmall = labelmatrix(nosmall)>0;
    % Get borders from label1c, then identify if any lie a putative split point 
    tmp = diagnos.label1c;
    tmp(~nosmall) = 0;
    tmp(tmp==0) = max(tmp(:))+1;
    border_mask = (tmp-imerode(tmp,ones(3)))>0;
    border_mask(diagnos.label1c==0) = 0;
    diagnos.edge_borders = (border_mask) + (diagnos.label1c>0);
    endpt_val = bwmorph(border_mask,'endpoints').*diagnos.cut_pts;
    endpt_val(imdilate(mask_cut,ones(5))) = 0;
    endpt_val(endpt_val==0) = nan;
    branch_pts = bwmorph(bwmorph(border_mask,'skel',inf),'branchpoints');
    cc_branch = bwconncomp(border_mask,8);
    avg_inflect = @(pix) nanmean(endpt_val(pix));
    num_branch = @(pix) sum(branch_pts(pix));
    num_end = @(pix) sum(~isnan(endpt_val(pix)));
    check_avg = cellfun(avg_inflect,cc_branch.PixelIdxList);
    check_branch = cellfun(num_branch,cc_branch.PixelIdxList);
    check_ep = cellfun(num_end,cc_branch.PixelIdxList);
    
    mask_tmp = diagnos.label1b>0;
    mask_tmp(cell2mat(cc_branch.PixelIdxList((check_avg>(p.NuclearInflection-30)) & (check_branch<1) & (check_ep==2))')) = 0;

    diagnos.label1b2 = bwlabel(mask_tmp,4);
    diagnos.label1c(diagnos.label1b2==0) = 0;
    pairs  = [diagnos.label1b2(:),diagnos.label1c(:)];
    [~,~,ic] = unique(pairs,'rows');
    diagnos.label1c = reshape(ic,size(diagnos.label1b))-1;
    cc_inflect = label2cc(diagnos.label1b2);

    % Bridge oversegmented nuclear subobjects (from edge-based divisions) together by shape
    diagnos.label1 = bridgenuclei(diagnos.label1c, cc_inflect, cutoff,p.ShapeDef, p.debug);
    
else % No objects were found- skip all these steps.
    diagnos.label1= diagnos.label1b;
end
%%

%%- - - - - - - - - - - - - - - - - - - Label2 - - - - - - - - - - - - - - - - - - - - - - -
% "Weak" objects missed by standard methods
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
    % Remove small objects from (merged) watershed
    areas1 = cell2mat(struct2cell(regionprops(diagnos.watershed_remainder,'Area')));
    diagnos.watershed_remainder(ismember(diagnos.watershed_remainder,find(areas1<cutoff.Area(1)))) = 0;
    diagnos.weak_ranked2 = rankpixels(diagnos.watershed_remainder, nucleus1); % Rerank in merged watershed
    diagnos.watershed_remainder((imdilate(diagnos.watershed_remainder,ones(3))-diagnos.watershed_remainder)>0) = 0;

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
    diagnos.label2 = bridgenuclei(diagnos.label2a,bwconncomp(diagnos.label2a>0,4),cutoff,p.ShapeDef, p.debug);
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
        vals = cat(2,4*ones(1,floor(0.08*length(locs))),3*ones(1,floor(0.12*length(locs))),2*ones(1,floor(0.15*length(locs))));
        vals = cat(2,vals,ones(1,length(locs)-length(vals)));
        ranked_image(locs(sort_order)) = vals;
    end
end


