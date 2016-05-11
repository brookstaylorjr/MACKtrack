function [ax_handle, fig_handle] = degradationPlot(timepoints,data, method)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% [ax_handle, fig_handle] = degradationPlot(timepoints,data, method)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% DEGRADATIONPLOT makes a bar chart showing relative contributions of synthesis and degradation at specific points
% in time
%
% timepoints   time of each measurement
% data         matrix of measurements (rows = timepoints, cols = no inhib | single inhib | dual inhib
% method       method of single inhibition ('syn' or 'deg' - default is 'deg')  
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% Check 3rd input
if nargin<3
    method = 'deg';
end
assert(max(strcmp(method,{'syn','deg'})>0),'Specified method must be either ''deg'' or ''syn''')
setcolors;

if strcmp(method,'deg')
    deg_x = data(:,1) - data(:,2);
    syn_x = data(:,2) - data(:,3);
    basal_x = data(:,3);   
else
    syn_x = data(:,2) - data(:,1);
    deg_x = data(:,3) - data(:,2);
    basal_x = data(:,3);
end

x_lim = [min(timepoints) - range(timepoints)*0.1 , max(timepoints) + range(timepoints)*0.1];
bar_width = abs(0.3*(timepoints(2)-timepoints(1)));


fig_handle = figure('PaperPositionMode','auto','Position',positionfig(300,150));
hold on
plot(x_lim, [0 0],':','color',colors.grays{4},'LineWidth',2); % Plot zero line
for i = 1:length(timepoints)
    x_val = [timepoints(i)-bar_width/2 timepoints(i)+bar_width/2 ...
        timepoints(i)+bar_width/2 timepoints(i)-bar_width/2];
    y_val1 =  [deg_x(i), deg_x(i), 0 0];
    y_val2 = [0 0 basal_x(i)+deg_x(i), basal_x(i)+deg_x(i)];
    y_val3 = [basal_x(i)+deg_x(i), basal_x(i)+deg_x(i),...
        basal_x(i)+syn_x(i)+deg_x(i), basal_x(i)+syn_x(i)+deg_x(i)];
    fill(x_val,y_val1,colors.trail{1},'EdgeColor',colors.trail{5},'LineWidth',1.5)
    fill(x_val,y_val2,colors.trail{5},'EdgeColor',colors.trail{5},'LineWidth',1.5)
    fill(x_val,y_val3,'w','EdgeColor',colors.trail{5},'LineWidth',1.5)
end
plot([0 1 2],x(ids,1),'--o','MarkerFaceColor','w','Color',colors.trail{5},'LineWidth',1.5)
hold off
set(gca,'XLim',[x_lim])
ax_handle = gca;