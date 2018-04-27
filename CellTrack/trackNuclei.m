function  [CellDataOut, queue_out] =  trackNuclei(queue_in,CellData,curr_frame, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [CellDataOut, queue_out] =  trackNuclei(queue_in,CellData,curr_frame, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% TRACKNUCLEI tracks nuclei across 2 frames (frame 1 to frame 2), using later frames for error checking.
%
% queue_in       stack of label matricies to be analyzed
% CellData       structure containing information/properties for all cells
% curr_frame     current image being analyzed (frame1 image number)
% p              parameters structure from setupTracking.m script
%
% CellDataOut    modified from CellData
% queue_out      queue_in, but with bottom frame relabeled with tracked objects
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%%

disp('Tracking decisions:')


% Define edge image for use later in decision-making
imgedge = false(size(queue_in(1).nuclei));
offset = p.ImageOffset{end}-p.ImageOffset{end-1};
row_low = ceil(max([1, 1+p.MinNucleusRadius+offset(1)]));
row_hi = floor(min([size(imgedge,1),size(imgedge,1)-p.MinNucleusRadius+offset(1)]));  
col_low = ceil(max([1, p.MinNucleusRadius+offset(1)]));
col_hi = floor(min([size(imgedge,2),size(imgedge,2)-p.MinNucleusRadius+offset(2)]));  
imgedge([1:row_low,row_hi:end],:) = 1; % Edge rows
imgedge(:,[1:col_low,col_hi:end]) = 1; % Edge cols

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
tmp.centroidx = vect_centroidx' - p.ImageOffset{end}(2);
tmp.centroidy = vect_centroidy' - p.ImageOffset{end}(1);
tmp.area = vect_area';
tmp.perimeter = vect_perimeter';    
tmp.obj(tmp.area==0) = 0;
if isfield(queue_in,'intensity') % Add in intensity data, if selected.
    tmp.intensity = queue_in(end).intensity(:);
end

labeldata = [CellData.labeldata, tmp]; % Add in new frame
CellData.blocks(CellData.FrameOut<curr_frame,:) = 0; % consistency check - make sure dead blocks are emptied.
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

% Sort links on distance travelled and similiarity (avg of perim/area change)
if ~isempty(links)
    [~,~,rnk1] = unique(links(:,5));
    [~,~,rnk2] = unique(links(:,6));
    [~,resolve_order] = sort((rnk1*2)+rnk2,'ascend');
    links = links(resolve_order,:);
end
% Rank/resolve links on distance travelled and similarity (average of perimeter/area changes)
while ~isempty(links)
    % Resolve link (unless it tries to link together two old blocks)
    if (~ismember(links(1,1),old_blocks(:,links(1,2)))) || (~ismember(links(1,3),old_blocks(:,links(1,4))))
        [links,blocks] = resolvelink(blocks,links,labeldata,p);
    else
        links(1,:) = [];
    end
end
%%
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
            disp(['Error: couldn''t match (old) block #',num2str(i),'[', num2str(CellData.blocks(i,:)),']'])
        end
        
    end
end
blocks(isnan(blocks(:,1)),:) = [];
blocks = blocks(:,2:end); % Drop the old frame


% if curr_frame==50
%     disp('stop')
% end

% Error checking, stage 1: look for paired false positive/false negative events
[blocks_pre, blocks, CellData] = robberhandling(blocks_pre, blocks, labeldata, CellData,p);


labeldata(1) = []; 

% Decision-making on PRE-EXISTING blocks
out_blocks = find(blocks_pre(:,1)==0);


while ~isempty(out_blocks)
    i = out_blocks(1);
    if blocks_pre(i,1) ~= 0
        out_blocks(1) = [];
        continue;
    end
    
    % - - - - - - - True negative - - - - - - - - - - - - - - - - - - - - - - 
    if sum(blocks_pre(i,:))==0 % If object doesn't ever reappear in stack again, let it die off
        % [Make sure we didn't kill the cell off already]
        if CellData.FrameOut(i) > curr_frame
            CellData.FrameOut(i) = curr_frame-1;
            disp(['Set cell #',num2str(i),'''s frame out: ',num2str(curr_frame-1)])
        end
        out_blocks(1) = [];
    else % False negative - get last known object position and interpolate
            disp(['False negative: adding #',num2str(i), ' into frame'])
            tmp_mask = false(size(queue_in(1).nuclei));
            frm0obj = CellData.blocks(i,1);
            frm1obj = length(labeldata(1).obj)+1;
            frmX = find(blocks_pre(i,:)>0,1,'first');
            frmXobj = blocks_pre(i,frmX);
            xpt1 = floor(CellData.labeldata(1).centroidx(frm0obj) + (labeldata(frmX).centroidx(frmXobj)...
                         - CellData.labeldata(1).centroidx(frm0obj))/frmX);
            ypt1 = floor(CellData.labeldata(1).centroidy(frm0obj) + (labeldata(frmX).centroidy(frmXobj)...
                        - CellData.labeldata(1).centroidy(frm0obj))/frmX);
            % Put image offset back in, cap to image size.
            y_place = ypt1 + p.ImageOffset{1}(1);
            x_place = xpt1 + p.ImageOffset{1}(2);
            x_place(x_place<1) = 1; x_place(x_place>size(queue_in(end).nuclei,2)) = size(queue_in(end).nuclei,2);
            y_place(y_place<1) = 1; y_place(y_place>size(queue_in(end).nuclei,1)) = size(queue_in(end).nuclei,1);
            rad1 = floor(sqrt(CellData.labeldata(1).area(frm0obj) + (labeldata(frmX).area(frmXobj)...
                - CellData.labeldata(1).area(frm0obj))/frmX)/pi);
            rad1(rad1<1) =1;
            tmp_mask(y_place,x_place) = 1;
            tmp_mask = imdilate(tmp_mask,diskstrel(rad1));
            queue_in(1).nuclei(tmp_mask) = frm1obj;
            % Add new object into block/labelinfo
            blocks_pre(i,1) = frm1obj;
            labeldata(1).obj = [labeldata(1).obj; frm1obj];
            labeldata(1).centroidx = [labeldata(1).centroidx; xpt1];
            labeldata(1).centroidy = [labeldata(1).centroidy; ypt1];
            labeldata(1).area = [labeldata(1).area; pi*rad1.^2];
            labeldata(1).perimeter = [labeldata(1).perimeter; 2*pi*rad1];
            labeldata(1).intensity = [labeldata(1).intensity; median(labeldata(1).intensity)];
            out_blocks(1) = [];
    end 
end
%%
% Decision-making on NEW blocks
blocks(blocks(:,1)==0,:) = []; % Don't classify blocks until cell shows in present frame
newblocks = [];
blockparent = [];
blockedge = [];
for i = 1:size(blocks,1)
    % PRELIMINARY CHECK: is new object in >66% of frames?
    if sum(blocks(i,:)>0) <= round(0.66*size(blocks,2))
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
        % Get all potential parents (within specified DriftDistance)
        p_obj = CellData.labeldata(1).obj(CellData.labeldata(1).obj>0);
        p_dist = sqrt( (CellData.labeldata(1).centroidx(p_obj) - labeldata(1).centroidx(blocks(i,1))).^2 + ...
            (CellData.labeldata(1).centroidy(p_obj) - labeldata(1).centroidy(blocks(i,1))).^2 );
        p_obj = p_obj(p_dist<=p.DriftDistance);
        p_idx = find(ismember(CellData.blocks(:,1),p_obj));

        % Ensure "child" isn't due to spurious false positive -> if obj in new block only coincides w/ obj in candidate "parent"
        % blocks for a couple frames, combine those blocks.
        sep = sum((double(repmat(blocks(i,:)>0,[length(p_idx), 1])) + double(blocks_pre(p_idx,:)>0))==1,2);
        if max(sep) > round(p.StackSize*0.66)
            ind1 = p_idx(find(sep==max(sep),1,'first'));
            row1 = blocks_pre(ind1,:);
            sub_row =  blocks(i,:);
            row1(row1==0) = sub_row(row1==0);
            blocks_pre(ind1,:) = row1;
            blocks(i,:) = nan;
            continue;
        end

        % Filter 1: parents must be at least 10 hrs old
        p_ages = ((curr_frame-CellData.FrameIn(p_idx))/p.FramesPerHour);
        p_ages(CellData.FrameIn(p_idx)<4) = 1000;
        p_ages(CellData.Parent(p_idx)==0) = 1000;
        if curr_frame<3
            p_ages(:) = 0;
        end
        p_idx(p_ages<10) = [];
        % Filter 2: sister/parent candidates cannot be short-lived
        p_idx(sum(blocks_pre(p_idx,:)==0,2)>=round(size(blocks_pre,2)*0.75)) = [];
        % Filter 3: sister + parent must be present (to compute "landing")
        p_idx(blocks_pre(p_idx,1)==0) = [];
        
        % Compute ideal "landing spot" for sister -> assume symmetrically-opposed division
        if ~isempty(p_idx)
            % (if applicable) Filter 4: eliminate any potential sisters that are v. different in intensity.
            if isfield(labeldata,'intensity')
                max_diff = prctile(abs(labeldata(1).intensity(blocks(i,1))-labeldata(1).intensity),50);
                elim = abs(labeldata(1).intensity(blocks(i,1))-labeldata(1).intensity(blocks_pre(p_idx,1))) > max_diff;
                if sum(elim)<numel(p_idx) % Make sure we don't eliminate all candidates...
                    p_idx(elim) = [];
                end
            end
            if length(p_idx)>1
                x_coord = 2*CellData.labeldata(1).centroidx(CellData.blocks(p_idx,1)) - labeldata(1).centroidx(blocks_pre(p_idx,1));
                y_coord = 2*CellData.labeldata(1).centroidy(CellData.blocks(p_idx,1)) - labeldata(1).centroidy(blocks_pre(p_idx,1));
                test1 = hypot(x_coord - labeldata(1).centroidx(blocks(i,1)),y_coord - labeldata(1).centroidy(blocks(i,1)));
                p_idx = p_idx(find(test1==min(test1),1,'first'));
            end

            % Create new lineages
            newblocks = cat(1,newblocks, blocks(i,:), blocks_pre(p_idx,:));
            blockedge = cat(1,blockedge, 0, 0);
            blockparent = cat(1,blockparent,p_idx, p_idx);
            blocks_pre(p_idx,:) = 0; % Zero parent
            CellData.blocks(p_idx,:) = 0; % Zero parent in orig data as well
            CellData.FrameOut(p_idx) = curr_frame-1; % Assign parent the proper frame out
            cellnum = size(blocks_pre,1)+length(blockedge)-1;           
            disp(['Created sisters (from #', num2str(p_idx),') '...
                '#',num2str(cellnum),' & #',num2str(cellnum+1),...
                ': [',num2str(newblocks(end-1,:)),'] and [',num2str(newblocks(end,:)),']'])
            
        % OPTION 3: if no suitable parent is found, we better be very sure this cell exists
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
CellDataOut.FrameOut = [CellData.FrameOut; (max(p.TimeRange))*ones(size(blockedge))];
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


function [blocks_pre, blocks, CellData] = robberhandling(blocks_pre, blocks, labeldata, CellData,p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% blocks_pre, blocks, labeldata, CellData] = robberhandling(blocks_pre, blocks, labeldata, CellData)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% ROBBERHANDLING looks for paired/nearby phenomena: a false positive, and a well-defined 
% "new" cell. If it seems likely that the mistakes are related, it will try to correct
% these trajectories.
%
% labeldata needs to have full information (i.e. prev frame hasn't been discarded yet)
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%%
fn_idx = find(blocks_pre(:,1)==0 & sum(blocks_pre,2)>0);
old_block = [CellData.blocks(:,1), blocks_pre];


new_block = [blocks_pre; blocks];
new_block = [zeros(size(new_block,1),1), new_block];


for i = fn_idx(:)'
    % Re-match our f.p. block to new blocks / preexisting blocks
    block_tmp = [old_block(i,1), zeros(1,size(blocks_pre,2))]; 
    allblock_tmp = new_block;
    allblock_tmp(i,:) = 0; % Keep block from re-matching to itself.
    links = linkblock(block_tmp, allblock_tmp, 1, labeldata, p);
    if isempty(links); continue; end
    [~,~,rnk1] = unique(links(:,5));
    [~,~,rnk2] = unique(links(:,6));
    [~,resolve_order] = sort((rnk1*2)+rnk2,'ascend');
    robber_idx = find(allblock_tmp(:,links(resolve_order(1),4))==links(resolve_order(1),3),1,'first');
    
     % Possibility 1: the f.p. matches to a new block - reassign for the next frame, then flag for resegmentation
     if robber_idx > size(blocks_pre,1)
        % (Ensure that matched object is extremely high quality!)
        if sum(new_block(robber_idx,:)==0)<3
            blocks_pre(i,1) = allblock_tmp(robber_idx,2);
            blocks_pre(i,2:end) = 0;
            blocks(robber_idx-size(blocks_pre,1),:) = 0;
            CellData.blocks(i,2:end) = 0;
            new_block(robber_idx,:) = 0;
            disp(['Reassigned F.N. cell #',num2str(i),' with a false positive obj that appeared'])            
        end
     % Possibility 2: the f.p. had its true trajectory "robbed" - see if robber matches to a new object.
     else
        robber_block = old_block(robber_idx,:);
        robber_block(2:end) = 0;
        allblock_tmp = new_block;
        allblock_tmp(robber_idx,:) = 0; % Keep block from re-matching to itself.
        links = linkblock(robber_block, allblock_tmp, 1, labeldata, p);
        if isempty(links); continue; end
        links(links(:,4)>2,:) = []; % Restrict this search to present frame
        if isempty(links); continue; end
        [~,~,rnk1] = unique(links(:,5));
        [~,~,rnk2] = unique(links(:,6));
        [~,resolve_order] = sort((rnk1*2)+rnk2,'ascend');
        newguy_idx = find(allblock_tmp(:,links(resolve_order(1),4))==links(resolve_order(1),3),1,'first');
        newguy_idx = newguy_idx - size(blocks_pre,1);

        % Step 3: perform resolution (if applicable)
        if newguy_idx > 0
            if sum(new_block(newguy_idx,:)==0)<3 % ensure new obj is H.Q.
            % All checks passed. Perform switcheroo - empty out remainder of block so it's rematched
            CellData.blocks(i,2:end) = [blocks_pre(robber_idx,1), zeros(1,p.StackSize-2)];
            CellData.blocks(robber_idx,2:end) = [blocks(newguy_idx,1), zeros(1,p.StackSize-2)];
            blocks_pre(i,:) = [blocks_pre(robber_idx,1), zeros(1,p.StackSize-1)];
            blocks_pre(robber_idx,:) = [blocks(newguy_idx,1), zeros(1,p.StackSize-1)];
            blocks(newguy_idx,:) = 0;
            disp(['Block ', num2str(robber_idx),' stole block ', num2str(i),'''s trajectory! Reverting.'])
            continue; % Skip the rest of the loop
            end
        end
     end
end


