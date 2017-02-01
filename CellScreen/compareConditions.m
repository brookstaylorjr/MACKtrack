%function [] = experimentSummary(AllData,measure_name,bin_lim)

if nargin<3
    pct = [0.01 92];
end
measure_name = 'MeanNuc1';
%%
% Collect data together
all_cond = fieldnames(AllData);
all_vals1 = [];
all_vals2 = {};
all_wells = {};
for i = 1:length(all_cond)
    if ~isfield(AllData.(all_cond{i}).Measurements,measure_name)
        error(['Error: Measurement field ''', measure_name,''' does not exist in one or more conditions.'])
    end
    all_vals1 = cat(1,all_vals1,AllData.(all_cond{i}).Measurements.(measure_name)(:));
    all_vals2= cat(1,all_vals2, {AllData.(all_cond{i}).Measurements.(measure_name)(:)});
    all_wells = 
end


%% Section 1: small-multiple histograms (by experiment)
bin_lim = prctile(all_vals1, pct);
bins = linspace(bin_lim(1),bin_lim(2),64);
n_cols = 6;
n_rows = ceil(length(all_vals2)/n_cols);

figure('Position',positionfig(1200,n_rows*150),'Name', ['Distributions of ''',measure_name,'''(all conditions)'])
ha = tight_subplot(n_rows,n_cols,[0.005 0.005]);
clr = linspecer(100);
all_y = zeros(size(ha));
for i = 1:length(ha)
    histogram(all_vals2{i},bins,'FaceColor',clr(98,:),'Parent',ha(i),'EdgeColor','none')
    set(ha(i),'XTIckLabel',{},'YTickLabel',{},'XGrid','on','XLim',bin_lim,'YGrid','on')
    if mod(i-1,n_cols) == 0
        ylabel(ha(i),measure_name,'FontSize',11,'FontWeight','normal','interpreter','none')
    end
    tmp = get(ha(i),'YLim'); 
    all_y(i) = tmp(2);    
end

for i = 1:length(ha)
    set(ha(i),'YLim',[0 max(all_y)])
    text(mean(bin_lim),max(all_y),[num2str(i),') ', all_cond{i},' ( n=',num2str(length(all_vals2{i})),' )'],...
        'HorizontalAlignment','center','VerticalAlignment','top','Parent',ha(i))
end


%% Section 2: bar plot (means, performed by CONDITION)





%% Section 3: bar plot (means, performed by CONDITION)



