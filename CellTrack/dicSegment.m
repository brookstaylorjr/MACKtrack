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
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% Setup Create basic masks/label matrices from phase/nuclear images
cell_mask = data.mask_cell>0;
nuc_mask = data.nuclei>0;

%% Inflection point splitting

% Find all candidate inflection points
diagnos.modifier_mask = findinflection(cell_mask,nuc_mask,p);
if nargout>1
    diagnos.all_inflection = diagnos.modifier_mask+nuc_mask+cell_mask;
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
%%
% Use straight edges to further break up mask
cell_mask = bwareaopen(cell_mask,round(pi*p.MinNucleusRadius.^2));
cell_broken = cell_mask&~imdilate(diagnos.edge_straight,ones(3));


obj1 = bwconncomp(cell_mask);
obj2 = bwlabel(bwareaopen(cell_broken,round(pi*p.MaxNucleusRadius.^2)));
obj2 = removemarked(obj2,nuc_mask,'keep');
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
keep_edges = imerode(keep_edges,ones(3));
diagnos.modifier_mask(keep_edges) = 1;

% Turn off pixels from modifier mask, do initial segmentation
cell_mask(diagnos.modifier_mask) = 0;
cell_mask = imopen(cell_mask,diskstrel(2));
output.img_straight = abs((image_in-prctile(image_in(:),0.02))/diff(prctile(image_in(:),[0.02 98])));
output.img_straight(output.img_straight<0) = 0; output.img_straight(output.img_straight>1) = 1;
output.img_straight(imdilate(diagnos.edge_straight,ones(3))) = 1;
output.img_straight(diagnos.edge_straight) = 0;

diagnos.image_segment = output.img_straight;
lambda = .02;
diagnos.seeds1 = IdentifySecPropagateSubfunction(double(data.nuclei),double(output.img_straight),cell_mask,lambda);

% Perform segmentation and correct
output.seeds2 = propagatesegment(diagnos.seeds1, data.mask_cell>0, image_in,...
    round(p.MinCellWidth/2), data.nuclei, lambda);

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

