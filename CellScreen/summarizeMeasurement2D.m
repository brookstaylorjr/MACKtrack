function [] = summarizeMeasurement2D(AllData, measure_x, measure_y, x_lim, y_lim)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  [] = summarizeMeasurement2D(measurement_struct)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SUMMARIZEMEASUREMENT2D shows a summary of a single measurement (e.g. 'MeanNuc1') across all conditons. Output graphs:
% -> scatterplot (dscatter) of all cells in a given condition
%
% INPUTS:
% AllData       master structure
% measure_x     names of measurement to put on  axis(e.g. 'MeanNuc1')
% measure_y     name of measurement to put on y axis (e.g. 'MeanNuc2')
% x_lim         graph x limits (percentiles - defaults to [0.1 97])
% y_lim         graph y limits (percentiles - defaults to [0.1 97])
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

if nargin<5
        y_lim = [0.01 97];
    if nargin<4
        x_lim = [0.01 97];
    end
end


% Collect data together (to set graph limits) 
all_cond = fieldnames(AllData);
all_vals_x = [];
all_vals_y = [];
for i = 1:length(all_cond)
    all_vals_x = cat(1,all_vals_x,AllData.(all_cond{i}).Measurements.(measure_x)(:));
    all_vals_y = cat(1,all_vals_y,AllData.(all_cond{i}).Measurements.(measure_y)(:));

end

[all_x] = restructuredata(AllData,measure_x);
[all_y] = restructuredata(AllData,measure_y);

%
x_rng = prctile(all_vals_x, x_lim);
y_rng = prctile(all_vals_y, y_lim);

n_cols = 6;
n_rows = ceil(length(all_x)/n_cols);

figure('Position',positionfig(1150,n_rows*200),'Name', ['''',measure_x,''' vs. ''', measure_y,''' (all conditions)'])
ha = tight_subplot(n_rows,n_cols,[0.01 0.01]);

colormaps=  loadcolormaps;
for i = 1:length(all_cond)
    drops = isnan(all_x{i}) | isnan(all_y{i});
    [~,h1] = dscatter2(all_x{i}(~drops),all_y{i}(~drops),'parent',ha(i));
    alpha(h1,0.4)
    r = corr(all_x{i}(~drops),all_y{i}(~drops));
    set(ha(i),'XLim',x_rng,'YLim',y_rng,'XTickLabel',{},'YTickLabel',{},'XGrid','on','Ygrid','on')
    text(mean(x_rng),max(y_rng),[num2str(i),') ', all_cond{i},' ( n=',num2str(length(all_x{i})),' )'],...
        'HorizontalAlignment','center','VerticalAlignment','top','Parent',ha(i),'BackgroundColor','w',...
        'Interpreter','none')
        text(max(x_rng),min(y_rng),['r = ',num2str(r)],...
        'HorizontalAlignment','right','VerticalAlignment','bottom','Parent',ha(i),'BackgroundColor','w',...
        'Interpreter','none')
end
colormap(colormaps.viridis(end:-1:1,:))
