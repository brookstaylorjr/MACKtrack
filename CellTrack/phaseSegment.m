function [label_out, data] = phaseSegment(queue, original, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% phaseSegment   Segment clusters of cells in the binary mask of a phase contrast image.
% (Based on CellProfiler's "propogate" algorithm, but improves result using contour points
% and intensity mask.)
%
%
% queue       tracking info structure: cell mask, nucleus label matrix, and supporting information
% original    original phase contrast image
% p           parameters structure from SetupTracking.m
%
% label_out   label matrix showing one cell per nucleus in nucLabel
% data        diagnostic information
%
% Subfunctions
% matchclosest.m, IdentifySecPropagateSubfunction.cpp (compiled w/ mex)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


%- - - - - - - - - - - - - - - - - - - SETUP - - - - - - - - - - - - - - - - - - - - - - -
% Create basic masks/label matrices from phase/nuclear images
SEsmall = strel('diamond',p.StrelSize); 
SEround = strel('disk',floor(p.MinCellWidth/2),4);
phase_mask = queue(1).Cell>0;

phase_mask2 = imopen(phase_mask,ones(5)); % Break very thin connections
phase_mask2 = imopen(phase_mask2,SEsmall);
cellClusters = bwconncomp(phase_mask2);
label_out = bwlabel(phase_mask2,4);
clusterBounds = regionprops(cellClusters,'BoundingBox');
halos = queue(1).halos;
nucLabel = queue(1).Nucleus;

% Cell-cluster skeleton
skeleton = bwmorph(phase_mask2,'skel','Inf');

% Mark skeleton point nearest to each nuclear centroid
nucCentroids = round(cell2mat(struct2cell(regionprops(nucLabel,'Centroid'))')');
nucCentroids = [nucCentroids(2,:);nucCentroids(1,:)];
skeletonMarkers = matchclosest(skeleton,nucCentroids,size(phase_mask));

% Setup images for diagnosis/ label matrix modification
data.marked_skeleton = zeros(size(phase_mask));
data.all_inflection = zeros(size(phase_mask));
data.modifier_mask = false(size(label_out));

%- - - - - - - - - - - - - - - - - - MAIN LOOP - - - - - - - - - - - - - - - - - - - - - -
for i= 1:cellClusters.NumObjects    
    try
        % If object contains more than one nuclei, pull it out as subset of image
        nucsWithin = unique(nucLabel(cellClusters.PixelIdxList{i}));
        nucsWithin(nucsWithin==0) = [];    
        numNucs = length(nucsWithin);
        
        if numNucs>1
            subset = clusterBounds(i).BoundingBox;  
            subsetRows = [max([floor(subset(2)),1]),  min([floor(subset(2))+ceil(subset(4)), size(phase_mask,1)]) ] ;
            subsetCols = [max([floor(subset(1)),1]) , min([floor(subset(1))+ceil(subset(3)), size(phase_mask,2)])];  
            
            % Get subset of other necessary data channels/masks 
            seed_reassign = label_out(subsetRows(1):subsetRows(2),subsetCols(1):subsetCols(2));
            seed_reassign(seed_reassign~=i) = 0;
            nucLabel2 = nucLabel(subsetRows(1):subsetRows(2),subsetCols(1):subsetCols(2));
            nucLabel2(seed_reassign~=i) = 0;
            skeleton2 = skeleton(subsetRows(1):subsetRows(2),subsetCols(1):subsetCols(2));
            skeleton2(seed_reassign~=i) = 0;
            skeletonMarkers2 = skeletonMarkers(subsetRows(1):subsetRows(2),subsetCols(1):subsetCols(2));
            skeletonMarkers2(seed_reassign~=i) = 0; % Remove markers from other objects
            
            if strcmp(p.ImageType,'phase')
                halos2 = halos(subsetRows(1):subsetRows(2),subsetCols(1):subsetCols(2));
            end
            
            % Check to see if object contains viable inflection points
            testSubset = seed_reassign>0;
            testStrel = diskstrel(ceil(p.MaxInflection/2));
            testOpened = bwareaopen(imopen(testSubset,testStrel),sum(sum(getnhood(testStrel))));
            if length(unique(bwlabel(testSubset))) == length(unique(bwlabel(testOpened)))
                continue
            end
            % Also drop objects that take up more than 25% of frame (high time cost)
            if (size(testSubset,1)*size(testSubset,2)) > ((p.ImageSize(1)*p.ImageSize(2))/4)
                continue
            end
            
            % Prune skeleton from ends until we "hit" the skeleton markers. (Get rid of misclassified endpoints.)
            endpoints = bwmorph(skeleton2,'endpoints');
            endpoints = endpoints&~bwareaopen(endpoints,2,4); % (make sure all endpoints = 1 pixel)
            numEndpoints = sum(endpoints(:));
            while numEndpoints > 0
                skeleton2(endpoints) = 0;
                endpoints = bwmorph(skeleton2,'endpoints');
                endpoints(~skeleton2) = 0;
                endpoints = bwmorph(endpoints,'shrink','Inf');
                endpoints(skeletonMarkers2) = 0;   
                numEndpoints = sum(endpoints(:));
            end  

            % Match each existing endpoint to a point on the perimeter, and drop those points (subdividing perimeter)
            perim = bwperim(seed_reassign>0);
            dropPoints = matchclosest(perim,bwmorph(skeleton2,'endpoints'),size(skeleton2));
            perim(imdilate(dropPoints,ones(3))) = 0;
            perimSegments = bwlabel(perim);
            
            % Get distance transform to each perimeter segment- sum the distance to closest two segments, and invert
            allDist = zeros([size(skeleton2),max(perimSegments(:))]);
            for currSeg = 1:max(perimSegments(:))
                allDist(:,:,currSeg) = bwdist(perimSegments==currSeg).*skeleton2;
            end
            [allDist,sortMat] = sort(allDist,3);
            dist2 = allDist(:,:,1) + allDist(:,:,2);
            dist2 = max(dist2(:))-dist2;
            dist2(~skeleton2) = 0;
            localMax = bwmorph(imregionalmax(dist2),'shrink','Inf');
            localMax(skeletonMarkers2) = 0;
            localMax(nucLabel2>0) = 0;
            branchPoints = false(size(skeleton2));
            
            
            % Make and fill in [n x n] connectivity matrix: we cycle through combinations and construct all direct paths between
            % connected nuclei
            if numNucs>2
                pathLengths = Inf.*ones(numNucs);
                savedPaths = cell(numNucs);
                % Separate skeleton into its "building blocks": nucPoints, branchPoints, and segments connecting them
                branchPoints = bwlabel(imdilate(bwmorph(skeleton2,'branchpoints'),ones(3)) & skeleton2);
                for branchPt = 1:max(branchPoints(:));
                    % Check to make sure that branch point contains at least 3 branches
                    if sum(sum(bwmorph(imdilate(branchPoints==branchPt,ones(3))&~(branchPoints==branchPt)&skeleton2,'shrink','Inf'))) < 3
                        branchPoints(branchPoints==branchPt)  = 0;
                    end
                end
                nucPoints = bwlabel(imdilate(skeletonMarkers2,ones(3))&skeleton2);
                segments = skeleton2;
                segments((nucPoints>0)|(branchPoints>0)) = 0;
                % Cycle through the object pairs: test if they are connected, and find minimal set of segments that connect them
                for obj1 = 1:numNucs-1
                    for obj2 = obj1+1:numNucs
                        skelPath = bwlabel(segments|(branchPoints>0)|(nucPoints==obj1)|(nucPoints==obj2));
                        a = skelPath(nucPoints==obj1);
                        b = skelPath(nucPoints==obj2);
                        if a(1)==b(1) % Nuclear objects are connected- find path.
                            skelPath(skelPath~=a(1)) = 0;
                            % Remove segments not in the connected object
                            newSegments = bwlabel(segments&(skelPath>0));
                            % Remove objects that have "naked" endpoints (i.e. don't have the nuclei in question on one end)
                            markers1 = bwmorph(skelPath>0,'endpoints')&~(nucPoints>0);
                            keepSegments = removemarked(newSegments,markers1);
                            skelPath((newSegments>0)&~(keepSegments>0)) = 0;
                            % Finally, remove any segments that do not interfere with connection
                            for seg1 = 1:max(keepSegments(:))
                                skelTest = skelPath>0;
                                skelTest(keepSegments==seg1) = 0;
                                skelTest = bwlabel(skelTest);
                                a = skelTest(nucPoints==obj1);
                                b = skelTest(nucPoints==obj2); 
                                 if a(1)==b(1) % Nuclear objects are still connected- eliminate segment
                                    skelPath(keepSegments==seg1) = 0;
                                 end
                            end
                        % Save path for later analysis; fill in info in pathLengths
                        finalPath = bwlabel(skelPath>0);
                        a =  finalPath(nucPoints==obj1);
                        pathLengths(obj1,obj2) = sum(sum(finalPath==a(1)));
                        savedPaths{obj1,obj2} = (finalPath==a(1));
                        end
                    end
                end

            else % Only 2 nuclei- no need to go through everything.
                savedPaths = cell(1);
                savedPaths{1,1} = skeleton2;
                pathLengths = 1;
            end
            
            %  Resolve singular inflection points in order of shortest path to longest path
            [sorted, resolveOrder] = sort(pathLengths(:),'ascend');
            resolveOrder(sorted==Inf) = [];
            offMask = false(size(skeleton2));
              
            for path1 = 1:length(resolveOrder)
                [r, c] = ind2sub(size(pathLengths),resolveOrder(path1));
                % Find the best candidate for inflection point: 
                inflection_candidates = (savedPaths{r,c}&localMax)&(~offMask);
                if sum(inflection_candidates(:))>0

                    [inflection_r,inflection_c] = find(inflection_candidates);
                    test_value = 0;
                    inflect_cut = [];

                    for pt = 1:length(inflection_r)
                        closestSegs = sortMat(inflection_r(pt),inflection_c(pt),1:2);
                        closestSegs = closestSegs(:);       
                        ind1 = matchclosest(perimSegments==closestSegs(1),[inflection_r(pt),inflection_c(pt)]');
                        ind2 = matchclosest(perimSegments==closestSegs(2),[inflection_r(pt),inflection_c(pt)]');
                        steps = max([abs(ind1(1)-ind2(1)), abs(ind1(2)-ind2(2))])+1;
                        test_cut = sub2ind(size(skeleton2),round(linspace(ind1(1),ind2(1),steps)),round(linspace(ind1(2),ind2(2),steps)));

                        if strcmp(p.ImageType,'phase')
                            % Make sure "cut" isn't larger than maximum allowed
                            if length(test_cut) > p.MaxInflection
                                halo_pct = 0;
                            else
                                halo_pct = sum(halos2(test_cut))/length(test_cut);
                            end

                            if (halo_pct > test_value) && (halo_pct>0.5)
                                inflect_cut = test_cut;
                                test_value = halo_pct;
                            end
                        else
                             inflect_cut = test_cut;
                        end
                    end

                    offMask(inflect_cut) = 1;
                end
            end

            % Update the label modifier_mask, as well as diagnostic images
            data.modifier_mask(subsetRows(1):subsetRows(2),subsetCols(1):subsetCols(2)) = offMask | data.modifier_mask(subsetRows(1):subsetRows(2),subsetCols(1):subsetCols(2));           
            marked_skeleton_subset = skeletonMarkers2+(2*skeleton2);
            marked_skeleton_subset((nucLabel2>0) & (marked_skeleton_subset==0)) = 1;
            marked_skeleton_subset(perimSegments>0) = 5;
            marked_skeleton_subset(branchPoints>0) = 6;
            data.marked_skeleton(subsetRows(1):subsetRows(2),subsetCols(1):subsetCols(2)) = marked_skeleton_subset + data.marked_skeleton(subsetRows(1):subsetRows(2),subsetCols(1):subsetCols(2));
            all_inflection_subset = dist2;
            all_inflection_subset(imdilate(localMax,ones(3))) = max(all_inflection_subset(:));
            data.all_inflection(subsetRows(1):subsetRows(2),subsetCols(1):subsetCols(2)) = data.all_inflection(subsetRows(1):subsetRows(2),subsetCols(1):subsetCols(2))+all_inflection_subset;
        end
    catch exception
        disp('Segmentation error:')
        disp(exception.message)
        disp(exception.stack(1))
        continue
    end
end


%- - - - - - - - - - - - - - - - - - PROPAGATE - - - - - - - - - - - - - - - - - - - - - -
% Turn off pixels from modifier mask, do initial segmentation
data.modifier_mask = removemarked(bwlabel(data.modifier_mask),nucLabel>0)> 0;
data.modifier_mask = imdilate(data.modifier_mask,ones(3));
phase_mask2(data.modifier_mask) = 0;

lambda = .02; % (Has very little effect on final result)

if strcmp(p.ImageType,'phase')
    % Use halos mask to possibly identify isolated cells (surrounded by halo)
    darkCells = (~halos)& phase_mask;
    darkCells = imopen(darkCells,SEround);
    darkLabeled = bwconncomp(darkCells);
    data.seed_label1 = nucLabel;
    for i = 1:darkLabeled.NumObjects
        a = unique(nucLabel(darkLabeled.PixelIdxList{i}));
        a(a==0) = [];
        if length(a)==1
            data.seed_label1(darkLabeled.PixelIdxList{i}) = a;
        end
    end
    data.seed_label2 = IdentifySecPropagateSubfunction(double(data.seed_label1),double(original),phase_mask2,lambda);
    data.seed_label2(isnan(data.seed_label2)) = 0;
else
       data.seed_label2 = IdentifySecPropagateSubfunction(double(nucLabel),double(original),phase_mask2,lambda);
end

% b) Correct for errors in round-cell segmentation.
borders = (imdilate(data.seed_label2,ones(3))-data.seed_label2)>0.5;
thin_pieces = data.seed_label2>0;
thin_pieces(borders) = 0;
thin_pieces = (thin_pieces)&~imopen(thin_pieces, SEround);
% Remove small noise from opened image
thin_pieces = bwareaopen(thin_pieces,4);
% Add back border that was subtracted earlier
thin_pieces = thin_pieces|(imdilate(thin_pieces,ones(3)) & borders & (data.seed_label2>0) );
% Now, reassign objects based on who they share longest border with.
label_dilated = data.seed_label2;
label_dilated(thin_pieces) = 0;
label_dilated = imdilate(label_dilated,ones(3));
reassignLabel = bwconncomp(thin_pieces);
data.seed_reassign = data.seed_label2;
for i = 1:reassignLabel.NumObjects
    newObj = label_dilated(reassignLabel.PixelIdxList{i});
    newObj(newObj==0) = [];
    newObj = mode(newObj);
    data.seed_reassign(reassignLabel.PixelIdxList{i}) = newObj;
end
data.seed_reassign(isnan(data.seed_reassign)) = 0;

% % Make sure no "islands" exist- all cells should be contiguous with their nuclei
data.seed_contig = data.seed_reassign;
data.seed_contig((imdilate(data.seed_contig,ones(5))-data.seed_contig)> 0) = 0; % Turns borders off
data.seed_contig = bwlabel(data.seed_contig>0,4);
data.seed_contig = removemarked(data.seed_contig,nucLabel>0,'keep');
data.seed_label3 = IdentifySecPropagateSubfunction(double(nucLabel),double(original),data.seed_contig>0,lambda);



% Do final assignment
label_out = IdentifySecPropagateSubfunction(double(data.seed_label3),double(original),phase_mask,lambda);
label_out(isnan(label_out)) = 0;