function [CellData, queue_out] = initializeCellData(queue_in,p,verbose)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% INITIALIZECELLDATA initializes CellData structure, creating initial set of "blocks" of nuclei to be tracked
%
% queue_in   n-deep stack of label matricies for cell nuclei
% p          parameters structure
% verbose      (optional) display link information
%
% CellData    structure with lineage information for each cell
% queue_out   modified queue_in, with bottom stack relabeled with initial set of cells
%
% Subfunctions
% linkblock.m
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if nargin<3
    verbose=0;
end

% Loop through queue: calculate regionprops and populate labeldata + initial blocks
labeldata = struct;
blocks = [];

for frm = 1:length(queue_in)
    % Get regionprops
    props = regionprops(queue_in(frm).nuclei,'Area', 'Centroid', 'Perimeter');
    tmpcell = struct2cell(props);
    vect_area = cell2mat(tmpcell(1,:));
    tmpmat = cell2mat(tmpcell(2,:));
    vect_centroidx = tmpmat(1:2:end);
    vect_centroidy = tmpmat(2:2:end);
    vect_perimeter = cell2mat(tmpcell(3,:));

    % Make 1st entry in labeldata
    labeldata(frm).obj = (1:length(props))';
    labeldata(frm).centroidx = vect_centroidx';
    labeldata(frm).centroidy = vect_centroidy';
    labeldata(frm).area = vect_area';
    labeldata(frm).perimeter = vect_perimeter';    
    labeldata(frm).obj(labeldata(frm).area==0) = 0;
    
    % Make starting set of blocks
    frmblocks = zeros(length(labeldata(frm).obj),length(queue_in));
    frmblocks(:,frm) = labeldata(frm).obj;
    frmblocks(sum(frmblocks,2)==0,:) = []; % Don't include all-zero blocks
    blocks = cat(1,blocks,frmblocks);
    
end

% Get initial set of links
links = [];
for i = 1:(size(blocks,1)-1)
    newlinks =linkblock(blocks(i,:), blocks, i+1, labeldata, p);
    links = cat(1,links,newlinks);
end

% Rank links on distance travelled and similarity (average of perimeter/area changes)
while ~isempty(links)
    [~,tmp] = sort(links(:,5),'ascend');
    [~,rnk1] = sort(tmp);
    [~,tmp] = sort(links(:,6),'ascend');
    [~,rnk2] = sort(tmp);
    [~,resolve_order] = sort((rnk1*2)+rnk2,'ascend');
    links = links(resolve_order,:);
    [links,blocks] = resolvelink(blocks,links,labeldata,p,verbose);   
end

% Of remainder, look for errors in list
for i = 1:size(blocks,1)
   if blocks(i,1) == 0 % Potential F.N.
       % NOTE: I'm not doing anything with this at the moment; I'd like to distinguish between new cells vs F.N. in frame 1.
       blocks(i,1) = NaN;
   else % Potential F.P.
        % Check1: is object in > 50% of frames?
        if (sum(blocks(i,:)>0)/size(blocks,2)) <= 0.5
            blocks(i,1) = NaN;
        end
   end
end
blocks(isnan(blocks(:,1)),:) = [];

% Wrap-up: combine blocks, assign CellData fields 
CellData.blocks = blocks;
CellData.labeldata = labeldata;
CellData.FrameIn = ones(size(CellData.blocks,1),1);
CellData.FrameOut = ones(size(CellData.blocks,1),1)*(max(p.TimeRange)-p.StackSize+1);
CellData.Parent = zeros(size(CellData.blocks,1),1);
CellData.Edge = false(size(CellData.blocks,1),1);

% Relabel bottom of stack to match existing block order, assign into queue_out
label_bottom = zeros(size(queue_in(1).nuclei));
for i = 1:size(CellData.blocks,1)
   label_bottom(queue_in(1).nuclei == CellData.blocks(i,1)) = i;
end
queue_out = queue_in;
queue_out(1).nuclei = label_bottom;