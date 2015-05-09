function [inflection] = findinflection(mask_cell, label_nuc, skeleton, markers, max_inflection, filter_size)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% FINDINFLECTION   Given a marked skeleton
%
% mask_cell        mask of cells (try to find inflection points here)
% label_nuc        label matrix of nuclei
% skeleton         pruned "backbone" of cells mask- find inflection points along this line
% markers          closest point to nuclei (on skeleton)
% max_inflection   maximum allowed length of inflection line
% filter_size      size of median filter
% 
%
% inflection   output structure with points and supporting info
%
% Subfuntions: matchclosest.m
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


% - - - - Prune skeleton from ends - - - -
endpoints = bwmorph(skeleton,'endpoints');
endpoints = endpoints&~bwareaopen(endpoints,2,4); % (make sure all endpoints = 1 pixel)
skeleton_old = false(size(skeleton));

% Stop when skeleton stabilizes
iter = 0;
while (~isequal(skeleton_old,skeleton)) && (iter<400)
    skeleton_old = skeleton;
    skeleton(endpoints) = 0;
    endpoints = bwmorph(skeleton,'endpoints');
    endpoints(markers) = 0;
    iter = iter+1;
end

inflection.pruned = skeleton;
inflection.perim = bwperim(mask_cell);

% - - - - Get distance transform to each perimeter segment- sum the distance to closest two segments, and invert - - - - 
dist1 = medfilt2(bwdist(inflection.perim),[filter_size,filter_size]).*inflection.pruned;
toofar = dist1>(max_inflection/2);
dist1 = max(dist1(:))-dist1;
dist1(~inflection.pruned) = 0;
localMax = bwmorph(imregionalmax(dist1),'shrink','Inf');
localMax(markers) = 0;
localMax(label_nuc>0) = 0;
localMax(toofar) = 0;
inflection.points = bwmorph(localMax,'shrink','Inf');