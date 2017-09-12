function [ha, plot_order] = ranksmult(graph_data, rankfactor, varargin)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% ha = ranksmult(graph_data, rankfactor)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% RANKSMULT will create a small-multiples line graph, where individual trajectories are
% ordered and sorted using a separately-passed value 
%
% INPUT (REQUIRED)
% graph_data    [n x m] matrix of n trajectories, each consisting of m dimensions/timepts
% rankfactor    [n x 1] vector, where each value corresponds to its respective trajectory
% xvect         (optional) time vector used in plotting
% INPUT (OPTIONAL)
% 'x'           x vector used in plotting (i.e. horizontal axis)
% 'YLim'        enforced Y limits of graph
%
% OUTPUT:
% ha            handles to tight_sublplot axes 
% plot_order    indicies of plotted individuals 
%-- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

%% Create input parser object, add required params from function input
p = inputParser;
% Required: graph_data
valid_data = @(x) assert((isnumeric(x)&&ismatrix(x)),...
    'graph_data must be a matrix (where rows correspond to a set measurements for one individual');
addRequired(p,'graph_data',valid_data);
% Required: rankfactor
valid_rank = @(x) assert((isvector(x)&&(length(x)==size(graph_data,1))),...
    'rankfactor must be vector (with one entry per individual in graph_data)');
addRequired(p,'rankfactor',valid_rank);



% Optional parameters
addParameter(p,'x',1:size(graph_data,2),@(x) assert(isvector(x)&&(length(x)==size(graph_data,2)),...
    'Length of X vector needs to match observations (columns) in graph_data'));
addParameter(p,'YLim',prctile(graph_data(:),[5 99]), @(x) assert(numel(x)==2,'YLim must be in form [y_min y_max]'));

% Parse parameters, assign to variables
parse(p,graph_data, rankfactor, varargin{:})
measure_bounds = p.Results.YLim;
xvect = p.Results.x;
%% Rank using rank factors; scale 0 to 100
[rank_val,idx] = sort(rankfactor,'ascend');

rank_val = floor((rank_val-min(rank_val))/(max(rank_val)-min(rank_val))*100);

% Set maximum number of subplots, then order individuals (even spacing) if over maximum
num_plots = 120;
if length(idx) <= num_plots
    plot_order = idx;
else
    plot_order = idx(floor(linspace(1,length(idx),num_plots)));
    rank_val = rank_val(floor(linspace(1,length(idx),num_plots)));
end


% Set graph characteristics
line_colors = cbrewer('div','Spectral',121);
line_colors(51:70,:) = []; % Drop light-colored rows
xpos = max(xvect)-0.02*(max(xvect)-min(xvect));
ypos =  max(measure_bounds) - 0.26*diff(measure_bounds);

fig = figure('name','smallmultiples');
set(fig,'Position',[500, 350, 876, 1000]);
ha = tight_subplot(15,8);


for i =1:length(plot_order)
    plot(ha(i),xvect,graph_data(plot_order(i),1:length(xvect)),...
        'Color',line_colors(rank_val(i)+1,:), 'LineWidth',2)
    set(ha(i),'XLim',[min(xvect) max(xvect)],'YLim',measure_bounds)
    set(ha(i),'XTickLabel',{[]}) 
    set(ha(i),'YTickLabel',{[]})
    disp_value = rankfactor(plot_order(i));
    if ((abs(disp_value) > 0.01) && (abs(disp_value) < 1000)) || (disp_value == 0)
        disp_str = num2str(round(disp_value*100)/100);
    else
        exp = floor(log10(abs(disp_value)));
        number = disp_value/(10^exp);
        disp_str = [num2str(round(number*100)/100),'e',num2str(exp)];
    end
    
    text(xpos,ypos,['#',num2str(i),': x = ',num2str(disp_str)],'Parent',ha(i),'HorizontalAlignment','right')
end
