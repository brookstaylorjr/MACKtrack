function [cc_out] = label2cc(label_in, compressed)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [cc_out] = label2cc(label_in, compressed)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% label2cc is the inverse of "labelmatrix"- it converts label_in to a bwconncomp structure
% (structure will be renumbered contiguously by default)
%
% INPUT:
% label_in       input label matrix
% compressed     boolean flag which sets conncomp as contiguous or not - defaults to true 
%
%
% OUTPUT:
% cc_out         output conncomp structure (identical to the output of 'bwconncomp')
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin < 2
    compressed = true;
end


% convert label_in to a conncomp struct
L = label_in(:);
idx = (1:numel(L)).'; % Linear indices as column vector
cc_out.PixelIdxList = accumarray(L(L~=0),idx(L~=0),[],@(x){x}); % Make the PixelIdxList
if compressed
    if ~all(ismember(1:max(L),L)) % If L wasn't monotonically increasing, make it that way.
        idx = unique(L(L~=0));
        cc_out.PixelIdxList = cc_out.PixelIdxList(idx);
    end
end
if ~iscell(cc_out.PixelIdxList)
    cc_out.PixelIdxList = {cc_out.PixelIdxList};
end

cc_out.NumObjects = length(cc_out.PixelIdxList);
cc_out.ImageSize = size(label_in);
cc_out.Connectivity = 4;
