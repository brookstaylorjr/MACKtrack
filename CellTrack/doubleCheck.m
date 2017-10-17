function [output, diagnos] = doubleCheck(data, p, queue_data)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [output, diagnos] = doubleCheck(data,cell_image,p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% DOUBLECHECK performs additional (application-specific) checking across the memory stack,
% before tracking/cell segmentation is run.
%
% INPUTS:
% data         structure with other output info (label_nuc and mask_cell)
% p            parameters struture
% queue_data   full stack of segmentations, etc - allows some rough history checking
%
% OUTPUTS:
% label_out    modified nuclear label
% data         output structure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% Output structure setup:
output.nuclei = data.label_nuc;
diagnos = struct;

% CHECK 1 (primary-specific): reassign cell masks to match identified nuclei
if strcmpi(p.ImageType,'none')
    % Overwrite cell masks
    output.mask0 = data.label_nuc>0;
    output.mask_cell = data.label_nuc>0;
end

% CHECK 2 (if selected): look for gross errors in nuclear segmentation (e.g. newest image is out of focus).
% Discard if really far out of range defined by prior stack.



% Save all information under diagnostic struct
diagnos = combinestructures(diagnos,output);

