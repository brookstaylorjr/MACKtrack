function  [CellDataOut, queue_out] =  trackNuclei(queue_in,CellData,curr_frame, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% TRACKNUCLEI tracks nuclei across 2 frames (frame 1 to frame 2), using later frames for error checking.
%
% queue_in       stack of label matricies to be analyzed
% CellData       structure containing information/properties for all cells
% currentImage   current image being analyzed (frame1 image number)
% p              parameters structure from setupTracking.m script
%
% CellDataOut    modified from CellData
% queue_out      queue_in, but with bottom frame relabeled with tracked objects
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
disp('Tracking decisions:')
% Define edge image for use later in decision-making
imgedge = false(size(queue_in(1).nuclei));
imgedge([1:p.MinNucleusRadius,end-p.MinNucleusRadius:end],:) = 1;
imgedge(:,[1:p.MinNucleusRadius,end-p.MinNucleusRadius:end]) = 1;

% Get regionprops of new frame (top of queue)
props = regionprops(queue_in(end).nuclei,'Area', 'Centroid', 'Perimeter');
tmpcell = struct2cell(props);
vect_area = cell2mat(tmpcell(1,:));
tmpmat = cell2mat(tmpcell(2,:));
vect_centroidx = tmpmat(1:2:end);
vect_centroidy = tmpmat(2:2:end);
vect_perimeter = cell2mat(tmpcell(3,:));

% Make new (temporary) structure, concatenate into old labeldata
tmp.obj = (1:length(props))';
tmp.centroidx = vect_centroidx';
tmp.centroidy = vect_centroidy';
tmp.area = vect_area';
tmp.perimeter = vect_perimeter';    
tmp.obj(tmp.area==0) = 0;
labeldata = [CellData.labeldata, tmp]; % Add in new frame
blocks = [CellData.blocks, zeros(size(CellData.blocks,1),1)]; % Add empty column to blocks

% Get rid of empty blocks
blocks(sum(blocks,2)==0,:) = [];
old_blocks = blocks;

% Get all unattached objects and add to current list of blocks
for frm = 2:length(queue_in)+1
    unattached = labeldata(frm).obj(~ismember(labeldata(frm).obj,blocks(:,frm)));
    unattached(unattached==0) = [];
    if ~isempty(unattached)
        new_blocks = zeros(size(unattached,1),size(blocks,2));
        new_blocks(:,frm) = unattached;
        blocks = [blocks; new_blocks];
    end
end

% Get initial set of links- don't link pre-existing objects together
links = [];
for i = 1:(size(blocks,1)-1)
    newlinks = linkblock(blocks(i,:), blocks, i+1, labeldata, p);
    links = cat(1,links,newlinks);
end

% Rank/resolve links on distance travelled and similarity (average of perimeter/area changes)
while ~isempty(links)
    [~,tmp] = sort(links(:,5),'ascend');
    [~,rnk1] = sort(tmp);
    [~,tmp] = sort(links(:,6),'ascend');
    [~,rnk2] = sort(tmp);
    [~,resolve_order] = sort((rnk1*2)+rnk2,'ascend');
    links = links(resolve_order,:);
    % Resolve link (unless it tries to link together two old blocks)
    if (~ismember(links(1,1),old_blocks(:,links(1,2)))) || (~ismember(links(1,3),old_blocks(:,links(1,4))))
        [links,blocks] = resolvelink(blocks,links,labeldata,p);
    else
        links(1,:) = [];
    end

end

% Pull pre-existing blocks (in numerically consistent order) 
blocks_pre = zeros(size(CellData.blocks));
for i = 1:size(CellData.blocks,1)
    col1 = find(CellData.blocks(i,:)>0,1,'first'); % Find 1st point where old block was defined
    if ~isempty(col1)
        obj = CellData.blocks(i,col1); % Get the object at that point 
        block1 = blocks(blocks(:,col1)==obj,2:end); % Find the (updated) block with that object 
        blocks(blocks(:,col1)==obj,:) = nan; % Move updated block from blocks to blocks_old
        try
            blocks_pre(i,:) = block1;
        catch ME
            disp(['Error in measurement:' , ME.message])
            disp('[',num2str(block1),']')   
            disp(['Error: couldn''t match (old) block #',num2str(i),'[', num2str(CellData.blocks(i,:)),']'])
        end
        
    end
end
blocks(isnan(blocks(:,1)),:) = [];
blocks = blocks(:,2:end); % Drop the old frame
labeldata(1) = []; 

% Decision-making on PRE-EXISTING blocks
mask_added = false(size(queue_in(1).nuclei));
for i = 1:size(blocks_pre,1)
    if blocks_pre(i,1) == 0
        if sum(blocks_pre(i,:))==0 % If object doesn't ever reappear in stack again, let it die off
            % Make sure we didn't kill the cell off already.
            if CellData.FrameOut(i) > curr_frame
                CellData.FrameOut(i) = curr_frame-1;
                disp(['Set cell #',num2str(i),'''s frame out: ',num2str(curr_frame-1)])
            end
        else % False negative; get last known object position and interpolate
            disp(['False negative: adding #',num2str(i), ' into frame'])
            tmp_mask = false(size(queue_in(1).nuclei));
            frm0obj = CellData.blocks(i,1);
            frm1obj = length(labeldata(1).obj)+1;
            frmX = find(blocks_pre(i,:)>0,1,'first');
            frmXobj = blocks_pre(i,frmX);
            xpt1 = floor(CellData.labeldata(1).centroidx(frm0obj) + (labeldata(frmX).centroidx(frmXobj) - CellData.labeldata(1).centroidx(frm0obj))/frmX);
            ypt1 = floor(CellData.labeldata(1).centroidy(frm0obj) + (labeldata(frmX).centroidy(frmXobj) - CellData.labeldata(1).centroidy(frm0obj))/frmX);
            rad1 = floor(sqrt(CellData.labeldata(1).area(frm0obj) + (labeldata(frmX).area(frmXobj) - CellData.labeldata(1).area(frm0obj))/frmX)/pi);
            tmp_mask(ypt1,xpt1) = 1;
            tmp_mask = imdilate(tmp_mask,diskstrel(rad1));
            mask_added = mask_added|tmp_mask; % Keep track of all added objects
            queue_in(1).nuclei(tmp_mask) = frm1obj;
            % Add new object into block/labelinfo
            blocks_pre(i,1) = frm1obj;
            labeldata(1).obj = [labeldata(1).obj; frm1obj];
            labeldata(1).centroidx = [labeldata(1).centroidx; xpt1];
            labeldata(1).centroidy = [labeldata(1).centroidy; ypt1];
            labeldata(1).area = [labeldata(1).area; pi*rad1.^2];
            labeldata(1).perimeter = [labeldata(1).perimeter; 2*pi*rad1];
        end 
    end
end

% Decision-making on NEW blocks
blocks(blocks(:,1)==0,:) = []; % Don't classify blocks until cell shows in present frame
newblocks = [];
blockparent = [];
blockedge = [];
for i = 1:size(blocks,1)
    % PRELIMINARY CHECK: is object in >60% of frames?
    if (sum(blocks(i,:)>0)/size(blocks,2)) <= 0.6
        blocks(i,:) = NaN;
    % OPTION 1: is object near the edge (drift-in)?
    elseif max(imgedge(queue_in(1).nuclei==blocks(i,1))) > 0 
        % Create new lineage
        newblocks = cat(1,newblocks,blocks(i,:));
        blockedge = cat(1,blockedge,1);
        blockparent = cat(1,blockparent,0);
        disp(['Drift-in cell #',num2str(size(blocks_pre,1)+length(blockedge)),': [',num2str(newblocks(end,:)),']'])
    
    % OPTION 2: division. Look at parent/sister candidates: can a suitable match be found?        
    else
        % Get 5 closest potential parents, then match/calculate distance to sisters
        candidates = CellData.labeldata(1).obj(CellData.labeldata(1).obj>0);
        parent_dist = sqrt( (CellData.labeldata(1).centroidx(candidates) - labeldata(1).centroidx(blocks(i,1))).^2 + ...
            (CellData.labeldata(1).centroidy(candidates) - labeldata(1).centroidy(blocks(i,1))).^2 );
        [~,p_sort1] = sort(parent_dist,'ascend');
        candidates = candidates(p_sort1(1:5));
        p_distances = parent_dist(p_sort1(1:5));
        blocknums = zeros(size(p_distances));
        for j = 1:length(blocknums)
            tmp = find(CellData.blocks(:,1)==candidates(j));
            if ~isempty(tmp)
                % Make sure (potential) sister is going to continue to exist
                if sum(blocks_pre(tmp,:)==0) < 2
                    blocknums(j) = tmp;
                end
            end
        end
        p_distances = p_distances(blocknums>0);
        blocknums = blocknums(blocknums>0);
        p_distances(blocks_pre(blocknums,1)==0) = [];
        blocknums(blocks_pre(blocknums,1)==0) = [];
        s_distances = sqrt( (labeldata(1).centroidx(blocks_pre(blocknums,1)) - labeldata(1).centroidx(blocks(i,1))).^2 + ...
            (labeldata(1).centroidy(blocks_pre(blocknums,1)) - labeldata(1).centroidy(blocks(i,1))).^2 );
        % For all candidates, parent should be CLOSER than sister
        blocknums = blocknums(p_distances<s_distances);
        p_distances = p_distances(p_distances<s_distances);
        % Of remaining candidates, minimize distance to parent
        p_block = blocknums(find(p_distances==min(p_distances),1,'first'));
        % a) Both daughters must exist in at least (n-1)/n frames 
        if ~isempty(p_block)
            test2a = max([sum(blocks_pre(p_block,:)==0), sum(blocks(i,:)==0)]) < 2;
        else
            test2a = 0;
        end
        % b) New daughter must be within 2 diameters of parent
        test2b =  min(p_distances) < (2*labeldata(1).perimeter(blocks(i,1))/pi);
        % c) If the parent has a parent, it must be at least 6 hrs old
        if CellData.Parent(p_block)==0
            test2c = 1;
        else 
            idx = CellData.Parent(p_block);
            test2c = ((curr_frame-CellData.FrameIn(idx))/p.FramesPerHour) > 6;
        end
        if test2a&&test2b&&test2c
            % Create new lineages
            newblocks = cat(1,newblocks, blocks(i,:), blocks_pre(p_block,:));
            blockedge = cat(1,blockedge, 0, 0);
            blockparent = cat(1,blockparent,p_block, p_block);
            blocks_pre(p_block,:) = 0; % Zero parent     
            CellData.FrameOut(p_block) = curr_frame-1; % Assign parent the proper frame out
            cellnum = size(blocks_pre,1)+length(blockedge)-1;
            disp(['Created sisters (from #', num2str(p_block),') '...
                '#',num2str(cellnum),' & #',num2str(cellnum+1),...
                ': [',num2str(newblocks(end-1,:)),'] and [',num2str(newblocks(end,:)),']'])
            
        % CHECK 4: if no suitable parent is found, we better be REAL sure this cell exists.
        elseif (sum(blocks(i,:)==0)) < 2
            % Create new lineage
            newblocks = cat(1,newblocks,blocks(i,:));
            blockedge = cat(1,blockedge,0);
            blockparent = cat(1,blockparent,0);
            cellnum = size(blocks_pre,1)+length(blockedge);
            disp(['Created drop-in cell #',num2str(cellnum),': [',num2str(newblocks(end,:)),']'])
        % Failing checks #2-4, drop the block.     
        else            
            blocks(i,:) = NaN;
        end
    end      
end

% Wrap-up: add in blocks and CellData fields
CellDataOut.FrameIn = [CellData.FrameIn; curr_frame*ones(size(blockedge))];
CellDataOut.FrameOut = [CellData.FrameOut; (max(p.TimeRange)-p.StackSize+1)*ones(size(blockedge))];
CellDataOut.Parent = [CellData.Parent; blockparent];
CellDataOut.Edge = [CellData.Edge; blockedge];
CellDataOut.labeldata = labeldata;
CellDataOut.blocks = [blocks_pre; newblocks];

% Relabel bottom of stack to match existing block order, assign into queue_out
label_bottom = zeros(size(queue_in(1).nuclei));
for i = 1:size(CellDataOut.blocks,1)
    if CellDataOut.blocks(i,1)>0
        label_bottom(queue_in(1).nuclei == CellDataOut.blocks(i,1)) = i;
    end
end
queue_out = queue_in;
queue_out(1).nuclei = label_bottom;
disp('- - - - - - - - - -')