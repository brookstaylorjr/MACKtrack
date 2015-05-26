function [CellData_out, queue_out] = memoryCheck(CellData, queue, cell_img, curr_frame, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% MEMORYCHECK checks each cell's history to identify inconsistencies in segmentation/tracking
% 
% CellData     flat structure with cell metadata and "blocks" used in future tracking
% queue        5-deep structure with cell and nuclear locations, as well as calculation intermediates
% cell_img     image of cells for re-segmentation
% curr_frame   current iteration of tracking
% p            tracking parameters structure
%
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
props_old= regionprops(queue(2).cells, 'Centroid');


% SWAP/JUMP check - make sure no cell leaps across another.
swap_cells = [];
swap_partners = [];
swap_displacement = [];
for n = reshape(checklist,1,length(checklist))
    % Get objects's X/Y pos in old/new frames
    x0=round(props_old(n).Centroid(1));
    y0=round(props_old(n).Centroid(2));
    x1=round(props_new(n).Centroid(1));
    y1=round(props_new(n).Centroid(2));
    % Generate its displacement line (convert to linear indicies)
    steps = max([abs(x0-x1), abs(y0-y1)])+1;
    displace_line = sub2ind(cc_all.ImageSize,round(linspace(y0,y1,steps)),round(linspace(x0,x1,steps)));
    displace_line(ismember(displace_line,cc_new.PixelIdxList{n})) = []; % Drop anything in same cell
    % See if most of displacement line runs through another cell    
    displace_line2 = displace_line(ismember(displace_line,cc_all.PixelIdxList{1}));
    if numel(displace_line2)/numel(displace_line) > 0.5
        swap_cells = cat(2,swap_cells,n);
        partner = mode(queue(1).cells(displace_line2));
        swap_partners = cat(2,swap_partners,partner);
        swap_displacement = cat(2,swap_displacement,sum(ismember(displace_line2,cc_new.PixelIdxList{partner})));
    end
end
% SWAP/JUMP resolution: classify error as either swap or jump, and fix appropriately
for i = 1:length(swap_cells)
    if swap_displacement(i) > (p.MinNucleusRadius*2)
        % Check to see if cell has matched partner in list, and if it also was a "big mover"
        test1 = swap_displacement(swap_cells==swap_partners(i)) > (p.MinNucleusRadius*2);
        test2 = swap_partners(swap_cells==swap_partners(i)) == swap_cells(i);
        test3 = ismember(swap_partners(i),checklist); % Partner needs have existed as well
        if isempty(test1) || isempty(test2)
            test1 = 0;
            test2 = 0;
        end
        if (length(test1)~=1) || (length(test2)~=1)
            test1 = 0;
            test2 = 0;
        end       
        if test1 && test2 && test3
            % Matched partner also moved. Swap positions and blocks.
            queue(1).cells(cc_new.PixelIdxList{swap_cells(i)}) = swap_partners(i);
            queue(1).cells(cc_new.PixelIdxList{swap_partners(i)}) = swap_cells(i);
            queue(1).nuclei(cc_nucs.PixelIdxList{swap_partners(i)}) = swap_cells(i);
            queue(1).nuclei(cc_nucs.PixelIdxList{swap_cells(i)}) = swap_partners(i);
            tmpblock = CellData.blocks(swap_cells(i),:);
            CellData.blocks(swap_cells(i),:) = CellData.blocks(swap_partners(i),:);
            CellData.blocks(swap_partners(i),:) = tmpblock;
            disp(['Switching positions of cells #', num2str(swap_cells(i)),' and #', num2str(swap_partners(i))]);
        else
            % No matched partner. Destroy old nucleus, flag it for re-adding, flag cell+partner for resegmentation
            queue(1).nuclei(cc_nucs.PixelIdxList{swap_cells(i)}) = 0;
            addlist = cat(2,addlist,swap_cells(i));
            fixlist = cat(2,fixlist,swap_cells(i),swap_partners(i));
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
            case 1 % Recreate dropped nucleus for resegmentation, edit block info
                r = mode(vals2(ismember(vals2,drops)));
                addlist = cat(2,addlist,r);
                fixlist = cat(2,fixlist,n,r);
                disp(['Cell # ',num2str(n),' error I - trying to add #',num2str(r)])   
            case 2 % Flag all implicated cells for resegmentation
                group1 = unique(vals2(~ismember(vals2, drops)));
                tmplist = [];
                for r = reshape(group1,1,length(group1));
                    if numel(vals2(vals2==r)) > (delta_a/4)
                        tmplist = cat(2,tmplist,r);
                    end
                end
                fixlist = cat(2,fixlist,[n,tmplist]);
                borderlist = cat(2,borderlist,[n,tmplist]);
                disp(['Cell # ',num2str(n),' error IIa - reassigning with cells [',num2str(tmplist),']'])
            case 3 % Remove the untracked cell from cell mask
                mask_reseg(untracked & queue(1).cells==n) = 0;
                fixlist = cat(2,fixlist,n);
                disp(['Cell # ',num2str(n),' error III - removing untracked cell from mask'])
            case 4 % Drop anything that was background in previous frame
                mask_reseg((queue(2).cells==0)&(queue(1).cells==n)) = 0;
                fixlist = cat(2,fixlist,n);
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
                r = mode(vals2(ismember(vals2,adds)));
                missing_frames = sum(sum(CellData.blocks([r,n],:)==0));
                if missing_frames > 0
                    % New nuc didn't pass; drop it and shift other label values down
                    droplist = cat(2,droplist,r);
                    disp(['Cell # ',num2str(n),' error V - cell #',num2str(r), ' FAIL (drop, shifting cells)'])
                    % If the dropped cell just "divided", also reassign its sister and drop it
                    parent = CellData.Parent(r);
                    if parent>0
                        x = find(CellData.Parent==parent,2,'last');
                        x(x==r) = [];   
                        droplist = cat(2,droplist,x);
                        addlist = cat(2,addlist,parent);
                        fixlist = cat(2,fixlist,parent);
                        disp(['   & "sister" cell (#', num2str(x),') reassigned back as parent #',num2str(parent)])
                    end
                    fixlist = cat(2,fixlist,n); % Add the kept cell to be fixed
                    CellData.blocks(n,2:end) = 0; % Drop remainder of existing cells's block
                else
                    disp(['Cell # ',num2str(n),' error V - cell #',num2str(r), ' PASS (keep)'])
                end
            case 2 % Fix segmentation error
                group1 = unique(vals2(~ismember(vals2, adds)));
                tmplist = [];
                for r = reshape(group1,1,length(group1));
                    if numel(vals2(vals2==r)) > (-delta_a/4)
                        tmplist = cat(2,tmplist,r);
                    end
                end
                fixlist = cat(2,fixlist,[n,tmplist]);
                borderlist = cat(2,borderlist,[n,tmplist]);
                disp(['Cell # ',num2str(n),' error IIb - reassigning with cells [',num2str(tmplist),']'])
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
        nuc_mask1 = (queue(2).nuclei==r) & (mask_reseg);
        if sum(nuc_mask1(:))>0 % Old nucleus must overlap with cells 
            nuc_mask = nuc_mask|nuc_mask1;
            nuc_seed(round(nprops_old(r).Centroid(2)), round(nprops_old(r).Centroid(1))) = r;
            % Keep blocks consistent: create new dummy obj
            new_ind = length(CellData.labeldata(1).obj)+1;
            CellData.labeldata(1).obj = [CellData.labeldata(1).obj; new_ind];
            CellData.blocks(r,1) = new_ind;
            CellData.FrameOut(r) = max(CellData.FrameOut);
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
            CellData.labeldata(1).centroidx = [CellData.labeldata(1).centroidx; new_props(r).Centroid(1)];
            CellData.labeldata(1).centroidy = [CellData.labeldata(1).centroidy; new_props(r).Centroid(2)];
            CellData.labeldata(1).area = [CellData.labeldata(1).area; new_props(r).Area];
            CellData.labeldata(1).perimeter = [CellData.labeldata(1).perimeter; new_props(r).Perimeter];
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
    mask_fix = ismember(queue(1).cells,fix1);
    mask_fix(all_cells>0) = 0;
    nucs_fix = queue(1).nuclei;
    nucs_fix(~ismember(nucs_fix,fix1)) = 0;
    nucs_old = label2cc(queue(2).nuclei,0);
    cells_old = label2cc(queue(2).cells,0);
    props_new = regionprops(nucs_fix,'Centroid');
    props_old = regionprops(nucs_old,'Centroid');
    cells_accum = zeros(size(mask_fix));
    try
        for i = 1:length(fix1)
            % Initialize cell mask; shift it by amount the nuclei putatively moved
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
    catch ME
        disp(['fix1 list: [',num2str(fix1(:)'),']']) 
        disp(getReport(ME,'extended'));
    end
    % And once more, making sure we got all those non-contiguous regions covered
    image_clamp = abs((cell_img-prctile(cell_img(:),0.02))/diff(prctile(cell_img(:),[0.02 98])));
    image_clamp(image_clamp<0) = 0; image_clamp(image_clamp>1) = 1;
    all_cells = IdentifySecPropagateSubfunction(all_cells,double(image_clamp),(queue(1).mask_cell),lambda);
    queue(1).cells = all_cells;
end

CellData_out = CellData;
queue_out = queue;



