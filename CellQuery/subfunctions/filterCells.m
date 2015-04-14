function DropCells = filterCells(CellData, FilterFlags, MinLifetime)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% DROPCELLS  Construct "DropCells", a logical array of cells that do not fit user-
% specified criteria

%         structure with user-specified criteria
% CellData        information passed from trackign process (see columns below) 
%
% FilterFlags:
% 1) non-edge 
% 2) parents
% 3) "doomed"
% 4) MinLifetime
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Columns in CellData: 
% 1) xyPosition
% 2) index in xy pos
% 3) FrameIn 
% 4) FrameOut 
% 5) Parent 
% 6) EdgeCell

dropInd = false(size(CellData(:,1)));

if FilterFlags(1) % Get vector of edge cells
    dropInd = dropInd | CellData(:,6);
end
if FilterFlags(2) || FilterFlags(3)% Get vector of non-parent cells and "doomed"cells
    parentsOnly = CellData(:,[1 5]);
    parentsOnly = parentsOnly(parentsOnly(:,2)>0,:);
    % Convert (e.g.) xy2 cell7 to a row index
    index_converter =  find(diff(CellData(:,1)));
    index_converter = cat(1,1,index_converter+1);
    parentRows = unique(index_converter(parentsOnly(:,1))+parentsOnly(:,2)-1);
    nonParents = true(size(CellData(:,1)));
    nonParents(parentRows) = 0;
    if FilterFlags(2)
        dropInd = dropInd | nonParents;
    end
    if  FilterFlags(3)
    % Define "doomed"cells as ones that did not divide, did not go off-edge, and whose end-frame is < max
    doomedCells = nonParents&(~CellData(:,6))&(CellData(:,4)~=0);
     dropInd = dropInd | ~doomedCells;
    end    
end
if FilterFlags(4) % Get vector of cells with lifetime under MinLifetime frames
    cellLifetimes = CellData(:,4)-CellData(:,3);
    cellLifetimes(cellLifetimes<0) = cellLifetimes(cellLifetimes<0)+max(CellData(:,4));
    dropInd = dropInd | cellLifetimes<MinLifetime;
end
DropCells = dropInd;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 