function [CellData_out, queue_out] = memoryCheck(CellData, queue, cell_img, curr_frame, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MEMORYCHECK checks each cell's history to identify inconsistencies in segmentation/tracking
% 
% CellData     flat structure with cell metadata and "blocks" used in future tracking
% queue        2-deep structure with past cell and nuclear locations, as well as calculation intermediates
% cell_img     image of cells for re-segmentation
% curr_frame   current iteration of tracking
% p            tracking parameters structure
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% ERROR TYPES:
%
% I.    A dropped cell was added to an existing cell; recreate the dropped cell
% IIa. Segmentation error with neighbors (area increase)
% IIb. Segmentation error with neighbors (area decrease)
% III. A perviously untracked cell was added
% IV.  Holes were improperly filled, increasing cell area
% V.   A false positive cell was added inappropriately
% VI.  Old filled holes were lost
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
disp('Memory checking:')
lambda = 0.02; 



% Data structures: bwconncomp, cell list
cc_old = label2cc(queue(2).cells,0);
cc_new = label2cc(queue(1).cells,0);
cc_all = label2cc(int8(queue(1).cells>0), 0);
cc_nucs = label2cc(queue(1).nuclei,0);
cells_new = unique(queue(1).cells(queue(1).cells>0));
cells_old = unique(queue(2).cells(queue(2).cells>0));
shift_old = queue(2).ImageOffset;
shift_new = queue(1).ImageOffset;
shift_tot = shift_new-shift_old;
image_size = cc_old.ImageSize;

% Data structures: masks/labels for resegmentation
mask_reseg = queue(1).cells>0;
fixlist = [];
borderlist = [];
addlist = [];
droplist = [];
% Get list of dropped cells and new cells
drops = find(CellData.FrameOut==(curr_frame-1));
adds = find(CellData.FrameIn==curr_frame);
checklist = cells_new(ismember(cells_new,cells_old));
% Calculate old and new centroids, and get conncomp of cell mask
props_new = regionprops(queue(1).cells, 'Centroid');
props_old= regionprops(queue(2).cells, 'Centroid','Area','Perimeter');



%% SWAP/JUMP check - make sure no cell leaps across another.
swap_cells = [];
swap_partners = [];
swap_displacement = [];
for n = reshape(checklist,1,length(checklist))
    if ~ismember(n,swap_partners)
        % Get objects's X/Y pos in old/new frames (in new frame's coordinates)
        [x0, y0] = getcentroid(n, props_old, shift_tot, image_size);
        [x1, y1] = getcentroid(n, props_new, [0 0], image_size);

        % Generate its displacement line (convert to linear indicies)
        steps = max([abs(x0-x1), abs(y0-y1)])+1;
        displace_line = sub2ind(cc_all.ImageSize,round(linspace(y0,y1,steps)),round(linspace(x0,x1,steps)));
        displace_line(ismember(displace_line,cc_new.PixelIdxList{n})) = []; % Drop anything in same cell
        % See if most of displacement line runs through another cell/nucleus    
        displace_line2 = displace_line(ismember(displace_line,cc_all.PixelIdxList{1}));
        if ~strcmp(p.ImageType,'none') && (numel(displace_line2)/numel(displace_line) > 0.5) % Cell images - use rel. measure
            jump_flag = 1;
        elseif strcmp(p.ImageType,'none') && (numel(displace_line2)>(1.5*p.MinNucleusRadius)) % Nuc images - use abs. measure
            jump_flag = 1;
        else
            jump_flag = 0;
        end
        if jump_flag
            swap_cells = cat(1,swap_cells,n);
            partner = mode(queue(1).cells(displace_line2));       
            swap_partners = cat(1,swap_partners(:),partner);
            swap_displacement = cat(1,swap_displacement(:),sum(ismember(displace_line2,cc_new.PixelIdxList{partner})));
        end
    end
end
%% SWAP/JUMP resolution: classify error as either swap or jump, and fix appropriately

% Reform old regionprops into structure
if ~isempty(swap_cells)
    tmpcell = struct2cell(props_old);
    vect_area = cell2mat(tmpcell(1,:));
    tmpmat = cell2mat(tmpcell(2,:));
    vect_centroidx = tmpmat(1:2:end);
    vect_centroidy = tmpmat(2:2:end);
    vect_perimeter = cell2mat(tmpcell(3,:));
    tmp.obj = (1:length(props_old))';
    tmp.centroidx = vect_centroidx' - queue(2).ImageOffset(2);
    tmp.centroidy = vect_centroidy' - queue(2).ImageOffset(1);
    tmp.area = vect_area';
    tmp.perimeter = vect_perimeter';    
    tmp.obj(tmp.area==0) = 0;
    labeldata_new = CellData.labeldata;
    if isfield(labeldata_new,'intensity'); labeldata_new = rmfield(labeldata_new,'intensity'); end

    % Classify each displacement as a JUMP or a SWAP.
    allblock = [zeros(size(CellData.blocks,1),1),CellData.blocks];
    swap_flag = zeros(size(swap_cells));
    for i = 1:length(swap_cells)
        % Check 1: was cell a big mover?
        if swap_displacement(i) < (p.MinNucleusRadius*2); continue; end

        % 2 flows - depends on if swap_partner is a new cell or not  
        if ismember(swap_partners(i),tmp.obj) % [Pre-existing cells]
            % Check 2: is swap_cell's highest preference (other than itself) its partner?
            block_tmp = [swap_cells(i), zeros(1,size(CellData.blocks,2))];
            allblock_tmp = allblock; allblock_tmp(swap_cells(i),:) = 0;
            links = linkblock(block_tmp, allblock_tmp, 1, [tmp,labeldata_new], p); 
            if isempty(links); continue; end
            [~,~,rnk1] = unique(links(:,5));
            [~,~,rnk2] = unique(links(:,6));
            [~,resolve_order] = sort((rnk1*2)+rnk2,'ascend');
            swap_idx = find(allblock(:,links(resolve_order(1),4))==links(resolve_order(1),3),1,'first'); 
            if swap_idx ~= swap_partners(i); continue; end
            % Check 3: is swap_partner's highest preference (other than itself) this cell?
            block_tmp = [swap_partners(i), zeros(1,size(CellData.blocks,2))];
            allblock_tmp = allblock; allblock_tmp(swap_partners(i),:) = 0;
            links = linkblock(block_tmp, allblock_tmp, 1, [tmp,labeldata_new], p); 
            if isempty(links); continue; end
            [~,~,rnk1] = unique(links(:,5));
            [~,~,rnk2] = unique(links(:,6));
            [~,resolve_order] = sort((rnk1*2)+rnk2,'ascend');
            swap_idx = find(allblock(:,links(resolve_order(1),4))==links(resolve_order(1),3),1,'first'); 
            if swap_idx ~= swap_cells(i); continue; end
            swap_flag(i) = 1;
        else
            % Check 2 (alt): is swap_cell's highest preference its partner?
            block_tmp = [swap_cells(i), zeros(1,size(CellData.blocks,2))];
            links = linkblock(block_tmp, allblock, 1, [tmp,labeldata_new], p); 
            if isempty(links); continue; end
            [~,~,rnk1] = unique(links(:,5));
            [~,~,rnk2] = unique(links(:,6));
            [~,resolve_order] = sort((rnk1*2)+rnk2,'ascend');
            swap_idx = find(allblock(:,links(resolve_order(1),4))==links(resolve_order(1),3),1,'first');
            if swap_idx ~= swap_partners(i); continue; end
            swap_flag(i) = 1;
        end
    end

    for i = 1:length(swap_cells)      
        if swap_flag(i)
            % A swap is more likely - switch trajectories of swapped cell and its partner
            queue(1).cells(cc_new.PixelIdxList{swap_cells(i)}) = swap_partners(i);
            queue(1).cells(cc_new.PixelIdxList{swap_partners(i)}) = swap_cells(i);
            queue(1).nuclei(cc_nucs.PixelIdxList{swap_partners(i)}) = swap_cells(i);
            queue(1).nuclei(cc_nucs.PixelIdxList{swap_cells(i)}) = swap_partners(i);
            tmpblock = CellData.blocks(swap_cells(i),:);
            CellData.blocks(swap_cells(i),:) = CellData.blocks(swap_partners(i),:);
            CellData.blocks(swap_partners(i),:) = tmpblock;
            disp(['Switching positions of cells #', num2str(swap_cells(i)),' and #', num2str(swap_partners(i))]);

        else % A jump was more likely. Destroy old nucleus, flag it for re-adding, flag cell+partner for resegmentation
            queue(1).nuclei(cc_nucs.PixelIdxList{swap_cells(i)}) = 0;
            addlist = cat(1,addlist(:),swap_cells(i));
            fixlist = cat(1,fixlist(:),swap_cells(i),swap_partners(i));
            checklist((checklist==swap_cells(i))|(checklist==swap_partners(i))) = [];
            % Finally, empty implicated blocks so that they're remade on next track cycle.
            CellData.blocks(swap_cells(i),2:end) = 0;
            CellData.blocks(swap_partners(i),2:end) = 0;
            disp(['Correcting cell #', num2str(swap_cells(i)),' jump over #',num2str(swap_partners(i))]);
        end
    end
end

% AREA CHANGE CHECK: look through remaining cells for large area increases/decreases
for n = reshape(checklist,1,length(checklist))        
    % Compare areas between consecutive frames
    delta_a =  (numel(cc_new.PixelIdxList{n}) - numel(cc_old.PixelIdxList{n}));
    pct_a = delta_a/numel(cc_old.PixelIdxList{n});
    
    % ERROR ID: 50%+ INCREASE
    if pct_a>(0.50)
        % Get source of new area
        vals = queue(2).cells(cc_new.PixelIdxList{n});
        vals2 = vals((vals~=0)&(vals~=n));
        % I. Dropped cells (nuclei ID error)
        error_nos(1) = numel(vals2(ismember(vals2,drops)));
        % II. Other cells (segmentation error)
        error_nos(2) = numel(vals2(~ismember(vals2,drops)));
        % III. Untracked cell (nuclei ID error)
        untracked = (queue(2).cells==0)&(queue(2).mask_cell);
        error_nos(3) = sum(untracked(cc_new.PixelIdxList{n}));
        % IV. Hole-filling error (cell ID error)
        new_pix = cc_new.PixelIdxList{n};
        new_pix(ismember(new_pix,cc_old.PixelIdxList{n})) = [];
        error_nos(4) = sum(~queue(1).mask0(new_pix));
        % Determine major source of increase
        [~,error_case] = max([(delta_a/2), error_nos]);
        error_case = error_case - 1;
        % ERROR HANDLING: INCREASES
        % disp([num2str(n) '''s error nos [',num2str(error_nos),'] - total error: ',num2str(delta_a)])
        switch error_case    
            case 1 % Area increase was from a previously dropped cell
                r = mode(vals2(ismember(vals2,drops)));
                % Dropped cell was parent - divvy up old parent to children 
                if ismember(r,CellData.Parent)
                    children = find(CellData.Parent==r);
                    children = children(:);
                    nucs = queue(1).nuclei;
                    nucs( (~ismember(nucs,children)) | (queue(2).cells~=r)) = 0;
                    if sum(ismember(children,unique(nucs(:)))) == length(children)
                        tmp = IdentifySecPropagateSubfunction(nucs,...
                            double(queue(1).img_straight),(queue(2).cells==r),lambda);
                        if sum(tmp(:)>0) == sum(queue(2).cells(:)==r)
                            queue(2).nuclei(queue(2).nuclei==r) = 0;
                            queue(2).nuclei = queue(2).nuclei + nucs;
                            queue(2).cells(queue(2).cells==r) = tmp(tmp>0);
                            fixlist = cat(1,fixlist(:),[n;children]);
                            borderlist = cat(1,borderlist(:),[n;children]); 
                            disp(['Cell # ',num2str(n),' error IIa ',...
                                '- reassigning with children [',num2str(children(:)'),']'])
                        end
                    end
                else 
                    % Recreate dropped nucleus for resegmentation, edit block info
                    addlist = cat(1,addlist(:),r);
                    fixlist = cat(1,fixlist,n(:),r);
                    disp(['Cell # ',num2str(n),' error I - trying to add #',num2str(r)])
                end
            case 2 % Flag all implicated cells for resegmentation
                group1 = unique(vals2(~ismember(vals2, drops)));
                tmplist = [];
                for r = reshape(group1,1,length(group1));
                    if numel(vals2(vals2==r)) > (delta_a/4)
                        tmplist = cat(1,tmplist(:),r);
                    end
                end
                fixlist = cat(1,fixlist(:),[n;tmplist]);
                borderlist = cat(1,borderlist(:),[n;tmplist]);
                disp(['Cell # ',num2str(n),' error IIa - reassigning with cells [',num2str(tmplist(:)'),']'])
            case 3 % Remove the untracked cell from cell mask
                mask_reseg(untracked & queue(1).cells==n) = 0;
                fixlist = cat(1,fixlist(:),n);
                disp(['Cell # ',num2str(n),' error III - removing untracked cell from mask'])
            case 4 % Drop anything that was background in previous frame
                mask_reseg((queue(2).cells==0)&(queue(1).cells==n)) = 0;
                fixlist = cat(1,fixlist(:),n);
                disp(['Cell # ',num2str(n),' error IV - removing filled holes from mask'])
        end
    end

    % ERROR ID: 25%+ DECREASE
    if pct_a<(-0.25)
        % Get source of lost area
        vals = queue(1).cells(cc_old.PixelIdxList{n});
        vals2 = vals((vals~=0)&(vals~=n));
        % V. New cell was added
        error_nos = [];
        error_nos(1) = -numel(vals2(ismember(vals2,adds)));
        % VI. Other cells (segmentation error)
        error_nos(2) = -numel(vals2(~ismember(vals2,adds)));
        % VII. Hole filling area
        pct1 = sum(~queue(1).mask0(cc_old.PixelIdxList{n})) / numel(cc_old.PixelIdxList{n}); 
        pct2 = sum(~queue(2).mask0(cc_new.PixelIdxList{n})) / numel(cc_new.PixelIdxList{n});
        error_nos(3) = delta_a*(pct2-pct1);
        % Determine major source of increase
        [~,error_case] = min([(delta_a/2), error_nos]);
        % disp([num2str(n) '''s error nos [',num2str(error_nos),'] - total error: ',num2str(delta_a)])
        error_case = error_case - 1;
        % ERROR HANDLING: DECREASE
        switch error_case
            case 1 % Double-check newly-added nucleus with stricter criteria
                tmp_subset = vals2(ismember(vals2,adds));
                r = mode(tmp_subset(:)); r = r(1);
                missing_frames = sum(sum(CellData.blocks([r,n],:)==0));
                if missing_frames > 0
                    % New nuc didn't pass; drop it and shift other label values down
                    droplist = cat(1,droplist,r);
                    disp(['Cell # ',num2str(n),' error V - cell #',num2str(r), ' FAIL (drop, shifting cells)'])
                    % If the dropped cell just "divided", also reassign its sister and drop it
                    parent = CellData.Parent(r);
                    if parent>0
                        x = find(CellData.Parent==parent,2,'last');
                        x(x==r) = []; % Drop self
                        droplist = cat(1,droplist(:),x(:));
                        addlist = cat(1,addlist(:),parent);
                        fixlist = cat(1,fixlist(:),parent);
                        disp(['   & "sister" cell (#', num2str(x(1)),') reassigned back as parent #',num2str(parent)])
                        if length(x)>1
                            disp(['NOTE: > 1 sister found: cells [',num2str(x(:)'),'] (??) - both were dropped.'])
                        end
                    end
                    fixlist = cat(1,fixlist,n); % Add the kept cell to be fixed
                    CellData.blocks(n,2:end) = 0; % Drop remainder of existing cells's block
                else
                    disp(['Cell # ',num2str(n),' error V - cell #',num2str(r), ' PASS (keep)'])
                end
            case 2 % Fix segmentation error
                group1 = unique(vals2(~ismember(vals2, adds)));
                tmplist = [];
                for r = reshape(group1,1,length(group1));
                    if numel(vals2(vals2==r)) > (-delta_a/4)
                        tmplist = cat(1,tmplist(:),r);
                    end
                end
                fixlist = cat(1,fixlist,[n; tmplist]);
                borderlist = cat(1,borderlist,[n; tmplist]);
                disp(['Cell # ',num2str(n),' error IIb - reassigning with cells [',num2str(tmplist(:)'),']'])
            case 3 % Add back anything from old label that is now called background
                queue(1).cells((queue(1).cells==0)&(queue(2).cells==n)) = n;
                disp(['Cell # ',num2str(n),' error VI - adding back lost area'])
        end
    end
end

% FIXING: curate lists
droplist = unique(droplist);
addlist = unique(addlist);

fixlist= unique(fixlist);
fixlist(ismember(fixlist,droplist)) = [];

borderlist = unique(borderlist);
borderlist(ismember(borderlist,droplist)) = [];
% FIXING: delete all cells that need to be deleted
for i = 1:length(droplist)
    r = droplist(i);
    CellData.blocks(r,:) = [];
    CellData.FrameIn(r) = [];
    CellData.FrameOut(r) = [];
    CellData.Parent(r) = [];
    CellData.Edge(r) = [];
    % Delete cell from mask, slide everything down
    queue(1).nuclei(queue(1).nuclei==r) = 0;
    queue(1).nuclei(queue(1).nuclei>r) = queue(1).nuclei(queue(1).nuclei>r)-1;
    queue(1).cells(queue(1).cells==r) = 0;
    queue(1).cells(queue(1).cells>r) = queue(1).cells(queue(1).cells>r)-1;
    % Also slide everything down in fixlist/borderlist
    fixlist(fixlist>r) = fixlist(fixlist>r)-1;
    borderlist(borderlist>r) = borderlist(borderlist>r)-1;
    addlist(addlist>r) = addlist(addlist>r)-1;
    droplist(droplist>r) = droplist(droplist>r)-1;
end

% FIXING: add back all cells that need to be added
if ~isempty(addlist) 
    % Initialize existing nuclei as points
    nprops_new = regionprops(queue(1).nuclei,'Centroid');
    nprops_old = regionprops(queue(2).nuclei,'Centroid');
    nucs = unique(queue(1).nuclei(queue(1).nuclei>0));
    nuc_seed = zeros(size(queue(1).nuclei));
    nuc_mask = queue(1).nuclei>0;
    for i = reshape(nucs,1,length(nucs))
        nuc_seed(round(nprops_new(i).Centroid(2)), round(nprops_new(i).Centroid(1))) = i;
    end
    % Loop addlist - make sure that nucleus is valid, add it as a new seed.
    for i = 1:length(addlist)
        r = addlist(i);     
        if ~strcmpi(p.ImageType,'none')
            nuc_mask1 = (queue(2).nuclei==r) & (mask_reseg);
        else
            nuc_mask1 = (queue(2).nuclei==r) | (mask_reseg); % need to make exception for nuclear-only segmentation
        end
        [ctr_c, ctr_r] = getcentroid(r,nprops_old, shift_tot, cc_old.ImageSize); % needs to be in cc_new's coords.
        
        if (sum(nuc_mask1(:))>0) && ~isnan(ctr_c) % Old nucleus must overlap with cells, and must be present in old frame
            nuc_mask = nuc_mask|nuc_mask1;           
            nuc_seed(ctr_r,ctr_c) = r;
            % Keep blocks consistent: create new dummy obj
            new_ind = length(CellData.labeldata(1).obj)+1;
            CellData.labeldata(1).obj = [CellData.labeldata(1).obj; new_ind];
            CellData.blocks(r,1) = new_ind;
            CellData.FrameOut(r) = max(CellData.FrameOut);
            disp(['Cell #', num2str(r),' passed - added back'])
        else
            disp(['No room to add add cell #', num2str(r), ' - deleting'])
            CellData.FrameOut(r) = curr_frame-1;
            CellData.blocks(r,:) = 0;
            borderlist(borderlist==r) = [];
            fixlist(fixlist==r) = [];
            queue(1).cells(queue(1).cells==r) = 0;
        end
    end
    % Get nuclei by using propagate function from seeds
    queue(1).nuclei = IdentifySecPropagateSubfunction(double(nuc_seed),double(nuc_mask),nuc_mask,0.02);
    % Now get region properties and reassign
    new_props = regionprops(queue(1).nuclei,'Area', 'Centroid', 'Perimeter');% Add nucleus props to labeldata
    for i = 1:length(addlist)
        r = addlist(i);
        if CellData.blocks(r,1)~=0
            CellData.labeldata(1).centroidx = [CellData.labeldata(1).centroidx; new_props(r).Centroid(1)-shift_new(2)];
            CellData.labeldata(1).centroidy = [CellData.labeldata(1).centroidy; new_props(r).Centroid(2)-shift_new(1)];
            CellData.labeldata(1).area = [CellData.labeldata(1).area; new_props(r).Area];
            CellData.labeldata(1).perimeter = [CellData.labeldata(1).perimeter; new_props(r).Perimeter];
            if isfield(CellData.labeldata,'intensity')
            CellData.labeldata(1).intensity = [CellData.labeldata(1).intensity; ...
                median(CellData.labeldata(1).intensity)]; % Fill in with dummy data
            end

        end
    end
end
% FIXING: resegment
if ~isempty(fixlist)
    all_cells = queue(1).cells;
    all_cells(queue(1).nuclei>0) = queue(1).nuclei(queue(1).nuclei>0);
    all_cells(ismember(all_cells,fixlist)) = 0;
    
    % Handle added/jumped cells
    borderlist(ismember(borderlist,addlist)) = [];
    fix1 = fixlist(~ismember(fixlist,borderlist));
    mask_fix = ismember(queue(1).cells,fix1);
    mask_fix(all_cells>0) = 0;
    nucs_fix = queue(1).nuclei;
    nucs_fix(~ismember(nucs_fix,fix1)) = 0;
    all_cells = all_cells + ...
        IdentifySecPropagateSubfunction(nucs_fix,double(queue(1).img_straight),mask_fix,lambda);
    % Handle mis-segmented cells: re-initialize old locations (in new mask)
    fix1 = borderlist;
    fix1(~ismember(fix1,queue(2).cells(:))) = [];
    mask_fix = ismember(queue(1).cells,fix1);
    mask_fix(all_cells>0) = 0;
    nucs_fix = queue(1).nuclei;
    nucs_fix(~ismember(nucs_fix,fix1)) = 0;
    nucs_old = label2cc(queue(2).nuclei,0);
    cells_old = label2cc(queue(2).cells,0);
    props_new = regionprops(nucs_fix,'Centroid');
    props_old = regionprops(nucs_old,'Centroid');
    cells_accum = zeros(size(mask_fix));
    for i = 1:length(fix1)
        try
            % Initialize cell mask for each cell; shift it by amount the nuclei putatively moved
            tmp = false(size(mask_fix));
            [r,c] = ind2sub(size(tmp),cells_old.PixelIdxList{fix1(i)});   
            r = round(r + props_new(fix1(i)).Centroid(2) - props_old(fix1(i)).Centroid(2));
            c = round(c + props_new(fix1(i)).Centroid(1) - props_old(fix1(i)).Centroid(1));
            r = max(r,1); r=min(r,size(tmp,1));
            c = max(c,1); c=min(c,size(tmp,2));
            tmp(sub2ind(size(tmp),r,c)) = 1;
            tmp = imerode(tmp,diskstrel(2));
            % Turn off any conflicting pixels, then add to mask.
            conflict = tmp&(cells_accum>0);
            cells_accum(conflict) = 0;
            tmp(conflict) = 0;
            cells_accum(tmp&mask_fix) = fix1(i);
        catch me
            disp(['cell #',num2str(fix1(i)),'triggered error:']) 
            disp(getReport(me,'extended','hyperlinks','off'));
            continue
        end
    end
    cells_accum(nucs_fix>0) = nucs_fix(nucs_fix>0);
    % Cycle through objects once; ensure every cell is contiguous
    cells_contig = zeros(size(cells_accum));
    for i = 1:length(fix1)
        tmp = cells_accum==fix1(i);
        locs = removemarked(bwconncomp(tmp,4),nucs_fix==fix1(i),'keep');
        cells_contig(cell2mat(locs.PixelIdxList')) = fix1(i);
    end
    cells_contig(all_cells>0) = 0;
    all_cells = all_cells + ...
        IdentifySecPropagateSubfunction(cells_contig,double(queue(1).img_straight),mask_fix,lambda);

    % Propagate outward once more; ensure all cell mask is accounted for, and that nuclei and cells match
    image_clamp = abs((cell_img-prctile(cell_img(:),0.02))/diff(prctile(cell_img(:),[0.02 98])));
    image_clamp(image_clamp<0) = 0; image_clamp(image_clamp>1) = 1;
    all_cells(~ismember(all_cells,unique(queue(1).nuclei(:)))) = 0;
    all_cells = IdentifySecPropagateSubfunction(all_cells,double(image_clamp),(queue(1).mask_cell),lambda);
    queue(1).cells = all_cells;
end

CellData_out = CellData;
queue_out = queue;


function [x_pos, y_pos] = getcentroid(n,rprops, shift_tot, image_size)

x_pos=round(rprops(n).Centroid(1))+shift_tot(2);
y_pos=round(rprops(n).Centroid(2))+shift_tot(1);
x_pos(x_pos<1)=1; x_pos(x_pos>image_size(2)) = image_size(2);
y_pos(y_pos<1)=1; y_pos(y_pos>image_size(1)) = image_size(1);


