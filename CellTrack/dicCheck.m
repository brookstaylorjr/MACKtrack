function [output, diagnos] = dicCheck(data,cell_image,p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [output, diagnos] = dicCheck(data,cell_image,p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% DICCHECK   Extract more information from original cell image to improve segmentation result
%
% INPUTS:
% data         structure with other output info (label_nuc and mask_cell)
% p            parameters struture
% aux          auxiliary info (TBD) 
%
% OUTPUTS:
% label_out    modified nuclear label
% data         output structure
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% - - - - Hough transform to find round cells in image - - - - 
houghSensitivity = 0.7;

% Gaussian-filter image
gauss_size = median(p.GaussSizes);
hsmall = gauss2D(gauss_size);
cell_image(cell_image==0) = min(cell_image(cell_image>0));
logimg = log(cell_image);
imfilt1 = imfilter(logimg, hsmall,'replicate');

% Hough- use size of largest holes/smallest cells. Sensitivity is fine on default.
radius1 = p.MinNucleusRadius;
if radius1>10
    [brightcenters, brightradii, brightfit] = imfindcircles(imfilt1,[radius1 radius1*3],'ObjectPolarity',...
        'dark','Sensitivity',houghSensitivity,'Method','twostage');
    diagnos.circleimg = zeros(size(imfilt1));
    output.circlelabel = zeros(size(imfilt1));
    for i = 1:length(brightfit)
        % Get circle center, dilate it (don't let it overlap over stronger circles)
        loc = [round(brightcenters(i,2)), round(brightcenters(i,1))];
        circle = diskstrel(floor(brightradii(i)));
        if (diagnos.circleimg(loc(1),loc(2)) == 0) && data.mask_cell(loc(1),loc(2))
            diagnos.circleimg(loc(1),loc(2)) = i;
            mask1 = imdilate(diagnos.circleimg==i,circle) & ((diagnos.circleimg==i)|(diagnos.circleimg==0));
            diagnos.circleimg(mask1) = brightfit(i);
            output.circlelabel(imerode(mask1,ones(3))) = i;
        end
    end


    % Fill in any strong (brightfit>0.15) circles that aren't already marked, ouput new nuclear label
    label_out = data.label_nuc;
    markers = bwareaopen((diagnos.circleimg>0.15)&data.label_nuc>0,p.NoiseSize);
    add_nucs = removemarked(output.circlelabel,markers,'remove');
    add_nucs(diagnos.circleimg<0.15) = 0;
    add_nucs(data.label_nuc>0) = 0;
    % Make sure new nuclei aren't super-close to edge of cell
    markers = ~imerode(data.mask_cell,ones(10))&(add_nucs>0);
    add_nucs = removemarked(add_nucs,markers,'remove');

    % Make a display figure
    diagnos.disp1 = cell_image;
    diagnos.disp1((imdilate(output.circlelabel>0,ones(3))-output.circlelabel) > 0) = min(cell_image(:));
    diagnos.disp1((imdilate(add_nucs,ones(3))-add_nucs) > 0) = max(cell_image(:));

    % Add the new nuclei to existing labelmatrix
    n = unique(add_nucs(add_nucs>0));
    for i = 1:length(n)
       label_out(add_nucs==n(i)) = max(data.label_nuc(:))+i;
    end
    output.nuclei = label_out;
else
    output.nuclei = data.label_nuc;
    diagnos = struct;
end
% Save all information under diagnostic struct
diagnos = combinestructures(diagnos,output);

