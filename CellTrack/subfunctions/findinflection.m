function [inflection] = findinflection(cell_mask, nuc_mask,p)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% [inflection] = findinflection(mask_cell, label_nuc, skeleton, markers, max_inflection)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% FINDINFLECTION   Given a marked skeleton
%
% INPUT:
% mask_cell        mask of cells (try to find inflection points here)
% label_nuc        label matrix of nuclei
% p                parameters structure from MACKtrack
% 
% OUTPUT:
% inflection      output matrix (sized as original image) with inflection lines marked
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Downsample image to a maximum size of 512x512
factor = round(max(size(cell_mask))/512);
cell_dwn = cell_mask(1:factor:end,1:factor:end);
nuc_dwn = nuc_mask(1:factor:end,1:factor:end);
max_inflection = round(p.MaxInflection/factor);
perim = bwperim(cell_dwn);

% Skeletonize cell clusters; eliminate larger blocks (2x3, 3x2, 3x3) from skeleton
skeleton = bwmorph(cell_dwn,'skel','Inf');
skeleton(perim) = 0;
check_ep = imopen(skeleton,ones(2)) & bwmorph(skeleton,'endpoints');
skeleton(bwmorph(check_ep,'shrink','inf')) = 0;

% Mark skeleton point closest to each nucleus
nucCentroids = round(cell2mat(struct2cell(regionprops(nuc_dwn,'Centroid'))'));
nucCentroids = [nucCentroids(:,2),nucCentroids(:,1)];
markers = matchclosest(skeleton,nucCentroids,size(cell_dwn));

% Prune skeleton from ends
endpoints = bwmorph(skeleton,'endpoints');
endpoints = endpoints&~bwareaopen(endpoints,2,4); % (make sure all endpoints = 1 pixel)
skeleton_old = false(size(skeleton));

iter = 0; % Stop when skeleton stabilizes
while (~isequal(skeleton_old,skeleton)) && (iter<400)
    skeleton_old = skeleton;
    skeleton(endpoints) = 0;
    endpoints = bwmorph(skeleton,'endpoints');
    endpoints(markers) = 0;
    iter = iter+1;
end

% Remove small holes in skeleton; save diagnostic outputs
skeleton = ~bwareaopen(~skeleton,4,4);
skeleton = bwmorph(skeleton,'skel','Inf');

% Get distance transform to nearest perimeter point; find regional minima
dist_img = bwdist(perim,'euclidean').*skeleton;
toofar = dist_img>(max_inflection/2);
dist_img = max(dist_img(:))-dist_img;
dist_img(~skeleton) = 0;
localMax = bwmorph(imregionalmax(dist_img),'shrink','Inf');
localMax(markers) = 0;
localMax(nuc_dwn>0) = 0;
localMax(toofar) = 0;

% Match each inflection point to its nearest perimeter point
[r,c] = find(localMax);
pts = [r(:), c(:)];
[r_perim,c_perim] = find(perim);
[matches, dists1] = matchclosest([r_perim(:),c_perim(:)],pts);

% Cycle matches; identify counterpart point on the other side of the skeleton, and match it to perimeter.
flipped_matches = pts + (pts-matches);
[matches2, dists2] = matchclosest([r_perim(:),c_perim(:)], flipped_matches);

% Checks: 
% - counterpart point is relatively close to perimeter
% - shape widens sufficiently around local neighborhood of point
nhood = 7;
ranges = zeros(size(pts,1),1);
dist_img(dist_img==0) = nan;
for i = 1:size(pts,1)
    subimg = dist_img(max([1,pts(i,1)-nhood]):min([size(dist_img,1),pts(i,1)+nhood]),...
        max([1,pts(i,2)-nhood]):min([size(dist_img,2),pts(i,2)+nhood]));
    ranges(i) = dist_img(pts(i,1),pts(i,2))-nanmin(subimg(:));
end
drop_pts = dists2 > max([2.24*ones(size(dists2)), 0.1*(dists1*2)],[],2);
drop_pts = (ranges<1.25)|drop_pts;
matches(drop_pts,:) = []; %dists1(drop_pts) = [];
matches2(drop_pts,:) = [];% dists2(drop_pts) = [];

%pts(drop_pts,:) = [];
%ranges(drop_pts) = [];
% In the future, we could improve this further - restrict inflection cuts to one per skeleton segment.
% (try to find point closest to middle of its two adjacent nucs?)


% Put matches into original image
matches_full = round(matches*factor);
matches2_full = round(matches2*factor);
perim_full = bwperim(cell_mask);
matches_full = matchclosest(perim_full,matches_full);
matches2_full = matchclosest(perim_full,matches2_full);


% Make inflection "cuts"
inflection = false(size(nuc_mask));
for pt = 1:size(matches2,1)
    steps = max([abs(matches_full(pt,1)-matches2_full(pt,1)), abs(matches_full(pt,2)-matches2_full(pt,2))])+1;
    test_cut = sub2ind(size(nuc_mask),round(linspace(matches_full(pt,1),matches2_full(pt,1),steps)),...
        round(linspace(matches_full(pt,2),matches2_full(pt,2),steps)));
    % Drop the point if it crosses nuclei
    if max(nuc_mask(test_cut))>0
        test_cut = [];
    end
    inflection(test_cut) = true;
end


% %%
% test = double(perim)-2;
% test(skeleton) = -1.2;
% test(sub2ind([512,512],matches(:,1),matches(:,2))) = -0.8;
% test(sub2ind([512,512],matches2(:,1),matches2(:,2))) = -0.4;
% for i = 1:length(ranges)
%     test(pts(i,1),pts(i,2)) = ranges(i);
% end
% test(imdilate(markers,ones(3))) = 0.2;
% figure,imagesc(test)



