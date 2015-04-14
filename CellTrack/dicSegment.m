function [output, diagnos] = dicSegment(data, original, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% DICSEGMENT   Segment clusters of cells in the binary mask of a phase contrast image.
% (Based on CellProfiler's "propogate" algorithm, but improves result using contour points)
%
% data       tracking info structure: cell mask, nucleus label matrix, and supporting information
% phaseOrig   original phase contrast image
% p           parameters structure from SetupTracking.m
%
% output      all information (masks) needed for tracking
% diagnos     major masks/thresholds created as intermediates
%
%
% Subfunctions
% matchclosest.m, inflectionID.m, IdentifySecPropagateSubfunction.cpp (compiled w/ mex)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

%% Setup Create basic masks/label matrices from phase/nuclear images
SEround = strel('disk',floor(p.MinCellWidth/2),4);
cell_mask = data.mask_cell>0;
cell_mask1 = imopen(cell_mask,SEround); % Break  thin connections
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
    inflection = findInflection(cell_dwn,nuc_dwn,skeleton,skeletonMarkers, max_inflection,p.MedianFilterSize);
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
    diagnos.disp_inflection = diagnos.modifier_mask+nuc_mask+cell_mask1;
catch exception
    disp('Segmentation error:')
    disp(exception.message)
    disp(exception.stack(1))
    disp(['(triggered by xy=',num2str(p.i),', t=',num2str(p.j),')'])
end
%% Propogation and correction

% Turn off pixels from modifier mask, do initial segmentation
diagnos.modifier_mask = removemarked(bwlabel(diagnos.modifier_mask),nuc_mask)> 0;
diagnos.modifier_mask = imdilate(diagnos.modifier_mask,ones(3));
cell_mask1(diagnos.modifier_mask) = 0;
lambda = .02; % (Has very little effect on final result)
diagnos.seeds1 = IdentifySecPropagateSubfunction(double(data.nuclei),double(original),cell_mask1,lambda);

% Perform segmentation and correct
output.seeds2 = propagateSegment(diagnos.seeds1, cell_mask, original, round(p.MinCellWidth/2), data.nuclei, lambda);

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

