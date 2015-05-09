function [curvature,idx] = calculatecurvature(object,L,n,sz)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

idx = sub2ind(sz,object(:,1),object(:,2));
idx = idx(1:length(unique(idx)));

object = object(1:length(idx),:);
curvature = nan(size(object,1),1);
for i = (n+1):(length(object)-n)
    delta_x = linspace(object(i-n,1),(object(i+n,1)),L) - object(i-n:i+n,1)';
    delta_y = linspace(object(i-n,2),(object(i+n,2)),L) - object(i-n:i+n,2)';
    curvature(i) = 1+sum((delta_x.^2+delta_y.^2).^0.5);
    %curv_disp(a(i,1),a(i,2)) = 1+sum((delta_x.^2+delta_y.^2).^0.5);
end
