function [output, diagnos] = doubleCheck(data, images, p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [output, diagnos] = doubleCheck(data,cell_image,p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% DOUBLECHECK performs additional (application-specific) checking across the memory stack,
% before tracking/cell segmentation is run.
%
% INPUTS:
% data         structure with other output info (label_nuc and mask_cell)
% images       structure with nuclear/cell images (if present) - images.nuc and images.cell
% p            parameters struture
% 
% queue_data   full stack of segmentations, etc - allows some rough history checking
%
% OUTPUTS:
% data         output structure (with modified labelmats and added info)
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


% Measure nuclear intensity data, if using.
if p.UseIntensity 
    tmpcc = label2cc(data.label_nuc,0);
    measurefunc = @(pix) nansum(images.nuc(pix));
    output.intensity = cellfun(measurefunc,tmpcc.PixelIdxList,'UniformOutput',1);    
end


% POSSIBLE FUTURE CHECK: look for gross errors in nuclear segmentation (e.g. newest image is out of focus).
% Discard if really far out of range defined by prior stack.



% Save all information under diagnostic struct
diagnos = combinestructures(diagnos,output);

