function cut_lines = perimetersplit(mask1,p)

% Set concave angle threshold (should be roughly 180+45 degrees)
angle_thresh = 220;


s = (p.MinNucleusRadius-4):(p.MinNucleusRadius-1);
s(s<1) = [];
b = bwboundaries(mask1,4);
sum_line = @(a,b) sum([a,b],2);
shifts = mat2cell(repmat(s(1),[length(b),1]),ones(size(b)));
[all_angles, ref_angles] = cellfun(@perimeterangles,b,shifts,'UniformOutput',0);
sum_angles = all_angles;
sum_refs = ref_angles;
for i = 2:length(s)
    shifts = mat2cell(repmat(s(i),[length(b),1]),ones(size(b)));
    [all_angles, ref_angles] = cellfun(@perimeterangles,b,shifts,'UniformOutput',0);
    sum_angles = cellfun(sum_line,all_angles,sum_angles,'UniformOutput',0);
    sum_refs = cellfun(sum_line,sum_refs,ref_angles,'UniformOutput',0);

end
rescale_line = @(a) a/length(s);
sum_angles = cellfun(rescale_line,sum_angles,'UniformOutput',0);
sum_refs = cellfun(rescale_line,sum_refs,'UniformOutput',0);

cut_lines = zeros(size(mask1));
% Threshold perimeter angles (per object)
for n = 1:length(sum_angles)
    idx = find(sum_angles{n}>angle_thresh);
    new_idx = [];
    dist = 8;
    while min(diff(idx)) < dist
        tmp = idx(find(sum_angles{n}(idx)==max(sum_angles{n}(idx)),1,'first'));
        new_idx = cat(1,new_idx,tmp);
        idx(ismember(idx,tmp-dist:tmp+dist)) = [];
    end
    idx = sort([new_idx;idx]);

    if length(idx)>1
        obj_angles = sum_angles{n}(idx);
        pts = b{n}(idx,:);
        v_angles =  sum_refs{n}(idx)/pi*180;

        bisect = v_angles+obj_angles/2;
        bisect(bisect>360) = bisect(bisect>360)-360;
        all_angles = zeros(size(pts,1),size(pts,2));
        for i = 1:size(pts,1)
            test_angles = atan2(-pts(:,1)+pts(i,1),pts(:,2)-pts(i,2))/ pi*180;
            test_angles(test_angles<0) = test_angles(test_angles<0)+360;
            all_angles(:,i) = test_angles-bisect(i);   
        end
        all_angles(sub2ind(size(all_angles),1:size(pts,1),1:size(pts,1))) = nan;
        all_angles = (abs(all_angles)+abs(all_angles)')/2;
        all_angles(all_angles>120) = nan;
        while sum(~isnan(all_angles(:)))>1
            [n1,n2] = find(all_angles==nanmin(all_angles(:)),1,'first');
            d1 = 1+max([abs(b{n}(idx(n1),1)-b{n}(idx(n2),1)),abs(b{n}(idx(n1),2)-b{n}(idx(n2),2))]);
            r = round(linspace(b{n}(idx(n1),1),b{n}(idx(n2),1),d1));
            c = round(linspace(b{n}(idx(n1),2),b{n}(idx(n2),2),d1));
            cut_lines(sub2ind(size(cut_lines),r,c)) = 2;     
            all_angles([n1 n2],:) = nan;
            all_angles(:, [n1 n2]) = nan;
        end
    end
end


function [perim_angles, ref_angles] = perimeterangles(pts,shift)

xy_pts = [pts(:,2), -pts(:,1)]; % Convert [r c] to [x y]
xy_fwd = circshift(xy_pts,shift,1);
xy_rev = circshift(xy_pts,-shift,1);
v1 = (xy_fwd - xy_pts)./repmat(sqrt(sum((xy_pts - xy_fwd).^2 , 2)),[1,2]);
v2  = (xy_rev - xy_pts)./repmat(sqrt(sum((xy_pts - xy_rev).^2 , 2)),[1,2]);
a1 = atan2(v2(:,2),v2(:,1)); a1(a1<0) = 2*pi+a1(a1<0);
a2 = atan2(v1(:,2),v1(:,1)); a2(a2<0) = 2*pi+a2(a2<0);
perim_angles = a1-a2; 
perim_angles(perim_angles<0) = perim_angles(perim_angles<0)+2*pi;
perim_angles(perim_angles>(2*pi)) = perim_angles(perim_angles>(2*pi))-(2*pi);
perim_angles = perim_angles/pi*180;
ref_angles = a2;





