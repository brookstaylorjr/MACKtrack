function [measure_out,shift_xy] = alignTrajectories(measurement, celldata, align_win, max_shift)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% ALIGNTRAJECTORIES takes input measurement (e.g. NFkB) and XYs by minimizing Euclidean
% distances. Trajectories will be shortened by 2*max_shift.
%
% [measure_out,shift_xy] = alignTrajectories(measurement, celldata, max_shift)
%
% measurement    input NxM matrix (M cells, N timepoint measurements)
% celldata       Nx6 cell information matrix (output from MACKtrack)
% align_win      (opt) window to use in calculating pairwise distances
% max_shift      (opt) maximum allowable shift for a given XY pos - defualts to 3
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

if nargin<4
    max_shift = 3;
    if nargin<3
        align_win = 60;
    end
end

xypos = unique(celldata(:,1));
M = size(measurement,2)-2*max_shift;
if align_win>M
    align_win = M;
end
orig = measurement;

% Assign 10th percentile of each trajectory to NaNs in that trajectory
lows = prctile(measurement,10,2);
reassign_matrix = repmat(lows,1,size(measurement,2));
measurement(isnan(measurement)) = reassign_matrix(isnan(measurement));


% Compute pairwise Euclidean distances between cells (using differing shifts)
avg_dists = zeros(length(xypos),length(xypos),max_shift+1);

% Get unshifted distances
d_nfkb = squareform(pdist(measurement(:,1:align_win),'euclidean'));
for i = 1:length(xypos)-1
    for j = i+1:length(xypos)
        rows_i = celldata(:,1) == xypos(i);
        rows_j = celldata(:,1) == xypos(j);
        avg_dists(i,j,1) = nanmean(nanmean(d_nfkb(rows_i,rows_j)));
        avg_dists(j,i,1) = nanmean(nanmean(d_nfkb(rows_i,rows_j)));
    end
end

% Get shifted distances
for i = 1:length(xypos)
    xy = xypos(i);
    rows_i = celldata(:,1)==xy;
    for fwd = 1:max_shift
        tmp = measurement(:,1:align_win);
        tmp(rows_i,:) = measurement(rows_i,1+fwd:align_win+fwd);
        d_nfkb = squareform(pdist(tmp,'euclidean'));
        % Get avg distance for fwd shifts, between current xy and rest of xys
        j_vect = find(xypos'~=xy);
        for j = j_vect
            rows_j = celldata(:,1) == xypos(j);
            avg_dists(i,j,fwd+1) = nanmean(nanmean(d_nfkb(rows_i,rows_j)));
        end
    
    end
    
end

[sorted, idx] = sort(avg_dists,3,'ascend');
sorted  = sorted(:,:,1);
idx = idx(:,:,1)-1;
idx = idx.*tril(ones(size(idx)))+idx.*triu(-ones(size(idx)));


test = cat(3,triu(sorted),triu(sorted.'));
idx2 = cat(3,triu(idx),triu(idx.'));

[~,idx3] = sort(test,3,'ascend');
idx3 = triu(idx3(:,:,1),1);

shifts = idx2(:,:,1).*(idx3==1) + idx2(:,:,2).*(idx3==2);
shifts = shifts - shifts.';
shift_xy = sum(shifts,2)/size(shifts,1);
shift_xy = round(shift_xy-max(shift_xy));

measure_out = orig(:,1:M);
for i =1:length(shift_xy)
    rows_i = celldata(:,1) == xypos(i);
    measure_out(rows_i,:) = orig(rows_i,1-shift_xy(i):M-shift_xy(i));
end

