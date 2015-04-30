function [output, diagnos] = dicSegment(data, image_in, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [output, diagnos] = dicSegment(data, original, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% DICSEGMENT segments clusters of cells in the binary mask of a phase contrast image.
% (Based on CellProfiler's "propogate" algorithm, but improves result using contour points)
%
% INPUT:
% data       tracking info structure: cell mask, nucleus label matrix, and supporting information
% phaseOrig   original phase contrast image
% p           parameters structure from SetupTracking.m
%
% OUTPUT:
% output      all information (masks) needed for tracking
% diagnos     major masks/thresholds created as intermediates
%
%
% Subfunctions
% matchclosest.m, findinflection.m, IdentifySecPropagateSubfunction.cpp (compiled w/ mex)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% Setup Create basic masks/label matrices from phase/nuclear images
SEround = strel('disk',floor(p.MinCellWidth/2),4);
cell_mask = data.mask_cell>0;
cell_mask1 = imopen(cell_mask,SEround); % Break thin connections
nuc_mask = data.nuclei>0;


%% Inflection point splitting

% Downsample image to a maximum size of 512x512
factor = round(max(size(cell_mask1))/512);
cell_dwn = cell_mask1(1:factor:end,1:factor:end);
nuc_dwn = data.nuclei(1:factor:end,1:factor:end);


% Skeletonize cell clusters; mark skeleton point closest to each nucleus
skeleton = bwmorph(cell_dwn,'skel','Inf');
skeleton(imopen(skeleton,ones(3))) = 0;
nucCentroids = round(cell2mat(struct2cell(regionprops(nuc_dwn,'Centroid'))')');
nucCentroids = [nucCentroids(2,:);nucCentroids(1,:)];
skeletonMarkers = matchclosest(skeleton,nucCentroids,size(cell_dwn));

diagnos.modifier_mask = false(size(cell_dwn));
try
    % Find all candidate inflection points
    max_inflection = round(p.MaxInflection/factor);
    inflection = findinflection(cell_dwn,nuc_dwn,skeleton,skeletonMarkers, max_inflection,p.MedianFilterSize);
    matches = matchclosest(inflection.perim, inflection.points);
    [r, c] = find(inflection.points);
    
    % Match each inflection point to its nearest perimeter point
    seeds = [r(:)';c(:)'];
    matches2 = round(matches - 2.5*(matches - seeds));
    matches2(matches2<1) = 1;
    matches2(1,matches2(1,:)>size(cell_dwn,1))= size(cell_dwn,1); 
    matches2(2,matches2(2,:)>size(cell_dwn,2))= size(cell_dwn,2); 
    % Create masks for modification/display
    for pt = 1:size(matches,2)
        % Drop the point if it doesn't clear other side of cell cluster
        if cell_dwn(matches2(1,pt),matches2(2,pt))==0   
            steps = max([abs(matches(1,pt)-matches2(1,pt)), abs(matches(2,pt)-matches2(2,pt))])+1;
            test_cut = sub2ind(size(cell_dwn),round(linspace(matches(1,pt),matches2(1,pt),steps)),...
                round(linspace(matches(2,pt),matches2(2,pt),steps)));
            % Drop the point if it crosses nuclei
            if max(nuc_dwn(test_cut))>0
                test_cut = [];
            end
            diagnos.modifier_mask(test_cut) = true;
        end
    end
    
    diagnos.modifier_mask = imresize(diagnos.modifier_mask,size(nuc_mask));
    diagnos.all_inflection = diagnos.modifier_mask+nuc_mask+cell_mask1;
catch exception
    disp('Segmentation error:')
    disp(exception.message)
    disp(exception.stack(1))
    disp(['(triggered by xy=',num2str(p.i),', t=',num2str(p.j),')'])
end

%% Straight-edge finding

% Form base edge mask from the default Canny edge filter - break/filter small branches
small = round(sqrt(p.NoiseSize));
edges = edge(image_in,'canny');
edges = bwmorph(edges,'skel',Inf);
edges([1,end],:) = 0;edges(:,[1,end]) = 0;
edge_broken = edges &~imdilate(bwmorph(edges,'branchpoints'),ones(3));
edge_broken = bwareaopen(edge_broken,small); % 

% Prune remaining skeleton 1x, then remove any remaining branch points
edge_pruned = edge_broken | (edges&imdilate(bwmorph(edges,'branchpoints'),ones(3)));
edge_pruned(bwmorph(edge_pruned,'endpoints')) = 0;
edge_pruned(imdilate(bwmorph(edge_pruned,'branchpoints'),ones(3))) = 0;
edge_pruned = bwareaopen(edge_pruned,small*2);

tic
% Trace edges and calculate curvature
[traced, edge_label] = bwboundaries(edge_pruned,8,'noholes');
L = 7; % (L must be odd)
n = (L-1)/2;
[curv,curv_idx] = cellfun(@calculatecurvature,traced,num2cell(L*ones(size(traced))),num2cell(n*ones(size(traced))),...
    num2cell(repmat(size(edge_pruned),length(traced),1),2),'UniformOutput',0);
curv_image = zeros(size(edge_pruned));
curv_image(cell2mat(curv_idx)) = cell2mat(curv);

% Drop edges where median curvature is too high
median_curv = cellfun(@nanmedian,curv);
label_list = 1:length(median_curv);
label_list(median_curv<n)=[];
edge_label(ismember(edge_label,label_list)) = 0;
diagnos.edge_straight = edge_label>0;

% Filter edges further (by nuclear overlap)
diagnos.edge_straight = diagnos.edge_straight & (curv_image<=L);
diagnos.edge_straight = removemarked(bwlabel(diagnos.edge_straight),imerode(nuc_mask,diskstrel(2))) > 0;
diagnos.edge_straight = diagnos.edge_straight | (curv_image>L);
diagnos.edge_straight = bwareaopen(diagnos.edge_straight,small);


% Use straight edges to remove offending inflection points
diagnos.modifier_mask = imdilate(diagnos.modifier_mask,ones(3));
diagnos.modifier_mask = removemarked(bwlabel(diagnos.modifier_mask),diagnos.edge_straight)> 0;
diagnos.modifier_mask = removemarked(bwlabel(diagnos.modifier_mask),nuc_mask)> 0;

% Use straight edges to potentially further break up mask
mask1 = bwareaopen(cell_mask,round(pi*p.MinNucleusRadius.^2));
mask2 = mask1&~imdilate(diagnos.edge_straight,ones(3));
obj1 = bwconncomp(mask1);
obj2 = bwlabel(bwareaopen(mask2,round(pi*p.MinNucleusRadius.^2)));
edgeobj = label2cc(imdilate(bwlabel(diagnos.edge_straight),ones(5)));
edgeobj2 = label2cc(imdilate(bwlabel(diagnos.edge_straight),ones(3)));

find_uniques = @(locs) reshape(unique(obj2(locs)),1,numel(unique(obj2(locs))));
unique_obj = cellfun(find_uniques,obj1.PixelIdxList,'UniformOutput',0);
test_length = @(vals) numel(vals(vals>0))>1;
idx = cellfun(test_length,unique_obj);
broken_obj = cell2mat(unique_obj(idx));
obj2(~ismember(obj2,broken_obj(broken_obj>0))) = 0;

find_uniques = @(locs) reshape(unique(obj2(locs)),1,numel(unique(obj2(locs))));
unique_obj = cellfun(find_uniques,edgeobj.PixelIdxList,'UniformOutput',0);
idx = cellfun(test_length,unique_obj);

keep_edges = cell2mat(edgeobj2.PixelIdxList(idx));
diagnos.modifier_mask(keep_edges) = 1;

%% Propogation and correction

% Turn off pixels from modifier mask, do initial segmentation
cell_mask1(diagnos.modifier_mask) = 0;
image_clamp = abs((image_in-prctile(image_in(:),0.02))/diff(prctile(image_in(:),[0.02 98])));
image_clamp(image_clamp<0) = 0; image_clamp(image_clamp>1) = 1;
image_clamp(imdilate(diagnos.edge_straight,ones(3))) = 1;
image_clamp(diagnos.edge_straight) = 0;

diagnos.image_segment = image_clamp;
lambda = .02;
diagnos.seeds1 = IdentifySecPropagateSubfunction(double(data.nuclei),double(image_clamp),cell_mask1,lambda);

% Perform segmentation and correct
output.seeds2 = propagatesegment(diagnos.seeds1, cell_mask, image_in, round(p.MinCellWidth/2), data.nuclei, lambda);

% Fill all holes that are internal to a single object
output.cells = output.seeds2;
cell_holes = imdilate(imfill(output.cells>0,'holes') &~ (output.cells),ones(3));
hole_cc = bwconncomp(cell_holes);
for i = 1:hole_cc.NumObjects
    hole_vals = unique(output.cells(hole_cc.PixelIdxList{i}));
    hole_vals(hole_vals==0) = [];
    if numel(hole_vals) == 1
        output.cells(hole_cc.PixelIdxList{i}) = hole_vals;
    end
end

% Save all information under diagnostic struct
diagnos = combinestructures(diagnos,output);

