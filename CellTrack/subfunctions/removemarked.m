function [output_obj] = removemarked(input_obj, markerMask, type)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% REMOVEMARKED Remove "marked" areas from a label matrix. Any object in the markerMask 
% that shares pixels with an object in an input label matrix (labelMarker) will cause that object to be removed 
% in a new labelMatrix. Works with bwconncomp structure or conventional label matrix.
%
% input_obj       label matrix/bwconncomp structure that defines pixel groups that are to be removed
% markerMask      the set of markers used to find candidate groups in labelMarker
% type            "remove"/"keep" - keep everything that is not marked, or marked - defaults to "remove"
%
% output_obj      label matrix/bwconncomp structure with marked objects turned to zeros (same as background)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if nargin<3
    type = 'remove';
end

% Convert label matrix to a bwconncomp to speed things up 
if ~isstruct(input_obj)
   input_cc = label2cc(input_obj);
else
    input_cc = input_obj;
end


remove_obj = [];
for i = 1:input_cc.NumObjects     
    if strcmp(type,'keep')
        if sum(markerMask(input_cc.PixelIdxList{i}))==0
        remove_obj = cat(2,remove_obj,i);
        end
    else
        if sum(markerMask(input_cc.PixelIdxList{i}))>0
        remove_obj = cat(2,remove_obj,i);
        end
    end
end
output_obj = input_cc;
output_obj.NumObjects = output_obj.NumObjects-length(remove_obj);
output_obj.PixelIdxList(remove_obj) = [];
    
if ~isstruct(input_obj)
   output_obj = labelmatrix(output_obj);
end