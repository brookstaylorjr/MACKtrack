function [cut_lines, all_pts] = perimetersplit(mask1,p)

% Set concave angle threshold (should be roughly 180+45 degrees)
angle_thresh = 225;
angle_thresh2 = 210; % More lenient threshold (used after smoothing angles)

s = (0:2) + min([2,round(p.MinNucleusRadius/4)]);
b = bwboundaries(mask1,8);
sum_line = @(a,b) sum([a,b],2);
shifts = mat2cell(repmat(s(1),[length(b),1]),ones(size(b)));
[ang, ref_angles] = cellfun(@perimeterangles,b,shifts,'UniformOutput',0);
sum_angles = ang;
sum_refs = ref_angles;
for i = 2:length(s)
    shifts = mat2cell(repmat(s(i),[length(b),1]),ones(size(b)));
    [ang, ref_angles] = cellfun(@perimeterangles,b,shifts,'UniformOutput',0);
    sum_angles = cellfun(sum_line,ang,sum_angles,'UniformOutput',0);
    sum_refs = cellfun(sum_line,sum_refs,ref_angles,'UniformOutput',0);
end
rescale_line = @(a) a/length(s);
sum_angles = cellfun(rescale_line,sum_angles,'UniformOutput',0);
sum_refs = cellfun(rescale_line,sum_refs,'UniformOutput',0);

% Threshold perimeter angles (per object)
cut_lines = false(size(mask1));
all_pts = zeros(size(mask1));
smooth_func = @(x) [x(end)+x(1)+x(2); x(1:end-2)+x(2:end-1)+x(3:end); x(end-1)+x(end)+x(1)]/3; % Circular running avg.
for n = 1:length(sum_angles)
    idx = find((sum_angles{n}>=angle_thresh) | (smooth_func(sum_angles{n})>=angle_thresh2));
    new_idx = [];
    dist = round(p.MinNucleusRadius/2);
    while min(diff(idx)) < dist
        tmp = idx(find(sum_angles{n}(idx)==max(sum_angles{n}(idx)),1,'first'));
        new_idx = cat(1,new_idx,tmp);
        idx(ismember(idx,tmp-dist:tmp+dist)) = [];
    end
    idx = sort([new_idx;idx]);
    all_pts(sub2ind(size(all_pts),b{n}(:,1),b{n}(:,2))) = sum_angles{n};

    if length(idx)>1
        obj_angles = sum_angles{n}(idx);
        pts = b{n}(idx,:);
        v_angles =  sum_refs{n}(idx)/pi*180;

        bisect = v_angles+obj_angles/2;
        bisect(bisect>360) = bisect(bisect>360)-360;
        ang = zeros(size(pts,1),size(pts,2));
        for i = 1:size(pts,1)
            test_angles = atan2(-pts(:,1)+pts(i,1),pts(:,2)-pts(i,2))/ pi*180;
            test_angles(test_angles<0) = test_angles(test_angles<0)+360;
            ang(:,i) = test_angles-bisect(i);   
        end
        ang(sub2ind(size(ang),1:size(pts,1),1:size(pts,1))) = nan;
        ang = (abs(ang)+abs(ang)')/2;
        
        count = 1;
        while sum(~isnan(ang(:)))>1
        % 1st cut should be between adjacent points -> maximize perim distance/cut distance ratio.
        % Subsequent cuts: minimize bisector angle difference (cuts should "point at" one another)
        if count==1
            d_perim = [diff(idx);length(b{n})-idx(end)+idx(1)];
            d_cut = b{n}(idx,:)-circshift(b{n}(idx,:),-1);
            d_cut = hypot(d_cut(:,1),d_cut(:,2));       
            n1 = find(d_perim./d_cut == max( d_perim./d_cut),1,'first');
            n2 = n1+1; n2(n2>length(idx)) = n2(n2>length(idx))-length(idx);
        else
            [n1,n2] = find(ang==nanmin(ang(:)),1,'first');
        end
            d1 = 1+max([abs(b{n}(idx(n1),1)-b{n}(idx(n2),1)),abs(b{n}(idx(n1),2)-b{n}(idx(n2),2))]);
            r = round(linspace(b{n}(idx(n1),1),b{n}(idx(n2),1),d1));
            c = round(linspace(b{n}(idx(n1),2),b{n}(idx(n2),2),d1));
            rm_idx = sub2ind(size(cut_lines),r,c);
            count = count+1;
            if (max(cut_lines(rm_idx))==0) && (min(mask1(rm_idx))==1) 
                cut_lines(rm_idx) = 1;     
                ang([n1 n2],:) = nan;
                ang(:, [n1 n2]) = nan;
            else
                ang(n1,n2) = nan;
                ang(n2,n1) = nan;
            end
                
        end
    end
end

%%
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





