function [] = summarizeMeasurement(AllData, measure_name, pct)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  [] = summarizeMeasurement(measurement_struct)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SUMMARIZEMEASUREMENT shows a summary of a single measurement (e.g. 'MeanNuc1') across all conditons. Output graphs:
% -> histogram of all cells in a given condition
% -> bar plot (by condition)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


if nargin<3
    pct = [0.01 97];
end

% Collect data together
all_cond = fieldnames(AllData);
all_vals1 = [];
for i = 1:length(all_cond)
    if ~isfield(AllData.(all_cond{i}).Measurements,measure_name)
        error(['Error: Measurement field ''', measure_name,''' does not exist in one or more conditions.'])
    end
    all_vals1 = cat(1,all_vals1,AllData.(all_cond{i}).Measurements.(measure_name)(:));
end

[by_condition, by_well, by_image] = restructuredata(AllData,measure_name);

%% Section 1: small-multiple histograms (by EXPERIMENT)
bin_lim = prctile(all_vals1, pct);
bins = linspace(bin_lim(1),bin_lim(2),64);
n_cols = 6;
n_rows = ceil(length(by_condition)/n_cols);

figure('Position',positionfig(1200,n_rows*150),'Name', ['Distributions of ''',measure_name,'''(all conditions)'])
ha = tight_subplot(n_rows,n_cols,[0.005 0.005]);
clr = linspecer(100);
all_y = zeros(size(ha));
for i = 1:length(ha)
    histogram(by_condition{i},bins,'FaceColor',clr(3,:),'Parent',ha(i),'EdgeColor','none')
    set(ha(i),'XTIckLabel',{},'YTickLabel',{},'XGrid','on','XLim',bin_lim,'YGrid','on')
    if mod(i-1,n_cols) == 0
        ylabel(ha(i),measure_name,'FontSize',11,'FontWeight','normal','interpreter','none')
    end
    tmp = get(ha(i),'YLim'); 
    all_y(i) = tmp(2);    
end

for i = 1:length(ha)
    set(ha(i),'YLim',[0 max(all_y)])
    text(mean(bin_lim),max(all_y),[num2str(i),') ', all_cond{i},' ( n=',num2str(length(by_condition{i})),' )'],...
        'HorizontalAlignment','center','VerticalAlignment','top','Parent',ha(i),'Interpreter','none')
end


%% Section 2: bar plots (means, performed by WELL, and means, performed by IMAGE)
colors = setcolors;
clr2 = clr(round(linspace(1,100,length(by_well))),:);

% Prep figure
figure('Position',positionfig(900,400));
ha= tight_subplot(2,1);


all_means = zeros(size(by_well));
all_err =  zeros(size(by_well));
for i = 1:length(by_well)
    m = cellfun(@nanmean,by_well{i});
    all_means(i) = mean(m);
    all_err(i) = std(m)/sqrt(length(m));
end

hold(ha(1),'on')
for i = 1:length(all_means)
    bar(i,all_means(i),'FaceColor',clr2(i,:),'EdgeColor','none','Parent',ha(1))
end
set(ha(1),'XLim',[0 length(all_means)*1.4],'XTick',[],'Box','on')
h = terrorbar(1:length(all_means),all_means,all_err,0.5);
set(h,'Color',colors.grays{3},'LineWidth',1,'Parent',ha(1));

hold(ha(1),'off')
text(0.5*mean(get(ha(1),'XLim')),max(get(ha(1),'YLim')),' Avearage across WELLS (w/ standard error)',...
    'HorizontalAlignment','center','VerticalAlignment','top','Parent',ha(1),'Interpreter','none')

% 2nd plot: group by image
all_means = zeros(size(by_image));
all_err =  zeros(size(by_image));
for i = 1:length(by_image)
    m = cellfun(@nanmean,by_image{i});
    all_means(i) = mean(m);
    all_err(i) = std(m);
end


hold(ha(2),'on')
for i = 1:length(all_means)
    bar(i,all_means(i),'FaceColor',clr2(i,:),'EdgeColor','none','Parent',ha(2))
end
set(ha(2),'XLim',[0 length(all_means)*1.4],'XTick',[],'Box','on')
h = terrorbar(1:length(all_means),all_means,all_err,0.5);
set(h,'Color',colors.grays{3},'LineWidth',1,'Parent',ha(2));
hold(ha(2),'off')
text(0.5*mean(get(ha(2),'XLim')),max(get(ha(2),'YLim')),' Avearage across IMAGES (w/ standard deviation)',...
    'HorizontalAlignment','center','VerticalAlignment','top','Parent',ha(2),'Interpreter','none')

l_ax = legend(all_cond,'Position',[0.7933  0.0758 0.15 0.9],'Interpreter','none');



