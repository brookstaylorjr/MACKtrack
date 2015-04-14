function links = linkblock(block1, blocks, start_pt, labeldata, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% LINKBLOCK is a tracking subfunction- it takes single block, and finds all possible blocks
% it could be linked to
%
% block1       incomplete row of cells (across queue of frames) to be matched
% block1props   area, perimeter, and centroid information for block1
% blocks       all candidate (incomplete) blocks that block1 can possibly match with
% blockprops   area, perimeter, and centroid information for all blocks
% start_pt     row index (of blocks) to begin search
% labeldata    structure with centroid/area/perimeter information, to which blocks correspond
% p            parameters structure (need p.DriftDistance)
%
% links       list of candidate links- gives pair of blocks, the frames that they crossover, 
%               and distance and similarity between them
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
links = [];
% Find frame range of block complementary blocks (i.e. ones that are zeroes where current block is defined)
comp_blocks = sum(blocks(:,block1~=0),2)==0;
if start_pt>1
    comp_blocks(1:(start_pt-1),:) = 0;
end

% Make sure we're not trying to match an empty block (or a full one)
if (sum(block1>0)==0) ||  (sum(block1>0)==length(block1))
    comp_blocks(:) = 0;
end

if max(comp_blocks) > 0
    % Find frames where block1 is defined
    loc_up = find(diff(block1>0)==1);
    loc_down = find(diff(block1>0)==-1);
    loc_up = [loc_up, length(block1)];

    % Get blocks who have frames in leading zeros of block1
    if isempty(loc_down) || (loc_up(1) < loc_down(1))
        ref_col = loc_up(1)+1;
        search_cols = 1:loc_up(1);
        search_blocks = (sum(blocks(:,search_cols),2)>0) & comp_blocks;
        newlinks = makelinks(block1, ref_col, search_cols(end:-1:1), blocks, search_blocks, labeldata,p);
        links = cat(1,links,newlinks);
        loc_up(1) = [];
    end
    % Cycle through loc_down and get other zero "subblocks" in
    for col = loc_down
        ref_col = col;
        search_cols = (col+1):min(loc_up(loc_up>col));
        search_blocks = (sum(blocks(:,search_cols),2)>0) & comp_blocks;  
        newlinks = makelinks(block1, ref_col, search_cols, blocks, search_blocks, labeldata,p);
        links = cat(1,links,newlinks);
    end
end
%===============================================================


function links = makelinks(block1, ref_col, col_vect, blocks, block_subset, labeldata, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MAKELINKS  given block1, its reference blocks
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
links = [];
for col = col_vect
    rows1 = block_subset & (blocks(:,col)>0);
    block_subset(rows1) = 0;
    dists = sqrt( (labeldata(col).centroidx(blocks(rows1,col)) - labeldata(ref_col).centroidx(block1(ref_col))).^2 + ...
        (labeldata(col).centroidy(blocks(rows1,col)) - labeldata(ref_col).centroidy(block1(ref_col))).^2 );
    % Turn off distances that are out of range
    rows1(rows1>0) = (dists<p.DriftDistance);    
    % Make links: [block1_obj, block1_col, block2_obj, block2_col, dist, delta_area/perim]
    newlinks = repmat([block1(ref_col), ref_col, 0, 0, 0, 0],sum(rows1),1);
    newlinks(:,3) = blocks(rows1,col);
    newlinks(:,4) = col;
    % Calculate distances and changes in area+perimeter
    newlinks(:,5) = sqrt( (labeldata(col).centroidx(blocks(rows1,col)) - labeldata(ref_col).centroidx(block1(ref_col))).^2 + ...
        (labeldata(col).centroidy(blocks(rows1,col)) - labeldata(ref_col).centroidy(block1(ref_col))).^2 );
    delta_area = (labeldata(col).area(blocks(rows1,col)) - labeldata(ref_col).area(block1(ref_col)) ) /...
        labeldata(ref_col).area(block1(ref_col));
    delta_perim = (labeldata(col).perimeter(blocks(rows1,col)) - labeldata(ref_col).perimeter(block1(ref_col)) ) /...
        labeldata(ref_col).perimeter(block1(ref_col));
    newlinks(:,6) = sqrt((delta_area.^2)+(delta_perim.^2));
    links = cat(1,links,newlinks);
end
