function [measure_out, divisions,divide_pts] = copychildren(measure_in, celldata_in, t)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [measure_out, divisions,divide_pts] = copychildren(measure_in, celldata_in)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% COPYCHILDREN uses lineage information to copy mother trajectories (i.e. prior to division) to any daughter cells.
% Requires a raw matrix of measurements (one row per cell), and "CellData" information matrix from MACKtrack.
% CellData info:  [ xy position | cell idx | frame in | frame out | parent | is_edge ]
%
% Optional argument: t, a time vector that corresponds to a measurement matrix (may include non-whole timepoints)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin<3
    t = 1:size(measure_in,2);
end

find_parent = @(row) find((celldata_in(:,1) == row(1)) & (celldata_in(:,2)== row(5)));
divisions = [];
divide_pts = false(size(measure_in));
measure_out = measure_in;
for i = 1:size(measure_out,1)
    if celldata_in(i,5)>0
        t_mom = (t>=1)&(t<=celldata_in(i,3));
        measure_out(i,t_mom) = measure_out(find_parent(celldata_in(i,:)),t_mom);
        divide_pts(find_parent(celldata_in(i,:)),celldata_in(i,3)) = 1;
        divide_pts(i,t_mom) = divide_pts(find_parent(celldata_in(i,:)),t_mom);
        divisions = cat(1,divisions,[i, celldata_in(i,3)]);
    end
end
