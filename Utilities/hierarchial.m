function [sortorder, links, p] = hierarchial(measurement, show_figure)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%HIERARCHIAL performs hierarchial clustering based on distance (not
% Euclidean) from 0
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Default is to hide figure
if nargin<2
show_figure = 0;
end


% Assign 10th percentile of each trajectory to NaNs in that trajectory
lows = prctile(measurement,10,2);
reassign_matrix = repmat(lows,1,size(measurement,2));
measurement(isnan(measurement)) = reassign_matrix(isnan(measurement));

% Do clustering/ordering based on Euclidean distance 
p  = pdist(measurement,'euclidean');
links = linkage(p); 

% Get sortorder from dendrogram (Supress initial dendrogram output)
tmpfig = figure('Visible','off');
C = cell(size(measurement,1),1);
C(:) = {''};
figure(tmpfig),[h,~,perm] = dendrogram(links,0,'orientation','left','labels',C);
if ~show_figure
close(tmpfig)
else
    dendrogram_axes = get(h(1),'Parent');
    set(dendrogram_axes,'LooseInset',get(dendrogram_axes,'TightInset')+ [0 0.03 0.02 0.01])
    set(h,'Color',[0.3 0.3 0.3])
    set(tmpfig,'Visible','on')
end
sortorder = perm(end:-1:1);
