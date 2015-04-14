function [trajectoriesOut, objects1,objects2,claim]  = overlapmatch(label1, label2,driftDist)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% OVERLAPMATCH:  Match "obvious" objects from one frame to the next using overlap method, then nearest neighbor.
%
% label1             first frame to be matched
% label2             second frame, objects map back to first
% driftDist          maximum allowable distance for one object to travel between 2 frames
%
% trajectoriesOut    (nx2) matrix (n is number of objects in frame 1) with matched objects
% objects1/objects2  remaining unmatched objects in frame1 and frame2 (respectively)
% claim              matrix with possible matches (in situations of non-unique/distinct overlap)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        

% - - - - - - 1. SETUP: Create lists of objects in frame1/ frame2, and trajectoriesOut matrix - - - - - - - - - - - - - - -
objects1 = unique(label1);
objects1(objects1==0) = [];
objects2 = unique(label2);
objects2(objects2==0) = [];

trajectoriesOut = zeros(length(objects1),2);
trajectoriesOut(:,1) = objects1;

% - - - - - - 2. OVERLAP: first (easy) categorization pass: look for distinct and unique, 1-to-1 overlap - - - - - - - - - - 
label1Locs = regionprops(label1,'PixelIdxList');
allMatches = [];
matches = cell(length(objects1),1);

% Cycle through label1 objects, look for overlaps
for j = 1:length(objects1)
    overlap = unique(label2(label1Locs(objects1(j)).PixelIdxList));
    overlap(overlap==0) = [];
    matches{j} = overlap;
    if ~isempty(overlap)
        allMatches= cat(2,allMatches,overlap');
    end
end
% Check to ensure each overlap is one-to-one
for j = 1:length(objects1)
    if (length(matches{j}) ==1) && (sum(allMatches==matches{j}) == 1)
        trajectoriesOut(j,2) = matches{j};
        objects1(objects1==trajectoriesOut(j,1)) = [];
        objects2(objects2==matches{j})=[];
    end
end

% - - - - - - 3). ORPHANS: match nearest remaining objects by distance between centroids - - - - - - - - - - - - - - -
orphans1 = zeros(size(label1));
orphans2 = zeros(size(label2));

if ~isempty(objects1)
    for objNo = 1:length(objects1)
        orphans1(label1== objects1(objNo)) =objects1(objNo);
    end
end
if ~isempty(objects2)
    for objNo =  1:length(objects2)
    orphans2(label2== objects2(objNo)) =objects2(objNo);   
    end
end
% Make claims matrix- fill with all possible claims and distance (precedence)
if ~isempty(objects1) && ~isempty(objects2)
    claim = []; 
    orphan1Props = regionprops(orphans1,'Centroid');
    orphan2Props = regionprops(orphans2,'Centroid');
    for i = 1:length(objects1)
      for j = 1:length(objects2)
          dist1 = norm(orphan1Props(objects1(i)).Centroid - orphan2Props(objects2(j)).Centroid);
          if dist1<driftDist
              claim = cat(1, claim, [objects1(i),objects2(j),dist1]);
          end    
      end
    end
    % Sort claims matrix
    if ~isempty(claim)
        [~,sortRows] = sort(claim(:,3));
        claim = claim(sortRows,:);
        % Go through claims matrix, resolving objects in order of lowest distance.
        for i = 1:size(claim,1)
                % Only match objects that have not been matched already
                if (ismember(claim(i,1),objects1) && ismember(claim(i,2),objects2))
                    trajectoryRow = trajectoriesOut(:,1) == claim(i,1);    
                    trajectoriesOut(trajectoryRow,2) = claim(i,2);
                    objects1(objects1 == claim(i,1)) = [];
                    objects2(objects2 == claim(i,2)) = [];
                end
        end
    end
else claim = []; 
end