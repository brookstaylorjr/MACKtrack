% Single-cell NFkB dynamics in a/b+e KO BMDMs

% Load single-cell datasets
if ~exist('nfkbKO.mat','file')
    nfkbKO = struct;
    ids = [                         
            292	% 2015-07-01_3.3ngTNF
            367	% 2015-12-08-aKO+3.3ngTNF
            371	% 2015-12-09_beKO+TNF (note: HETEROZYGOUS)
            ];
    for i = 1:length(ids)
        [nfkbKO(i).metrics, nfkbKO(i).fourier] = nfkbmetrics(ids(i));
    end
    save('nfkbKO.mat', 'nfkbKO')
else
    load('nfkbKO.mat')
end
savedir = '/Users/brooks/Dropbox/2016_R01_macrophage/figures/1.2/';
savestem = 'abeKO';
setcolors;loadcolormaps;

%%
num_cells = 200; % Number of rows to show
graph_lim = {[0.6 6] , [0.6 6], [0.3 3]};
t_max = 10;
osc_cutoff = 0.43; % Frequency cutoff for cell to be considered 'oscillatory'

% 3 x 1 heatmaps
figs.heatmaps(1) = figure('name','ColormapStack','Position', positionfig(520,210),'PaperPositionMode','auto');
ha = tight_subplot(1,3,[0.02 0.02]);
for i = 1:length(ha)
    x = medfilt1(nfkbKO(i).metrics.time_series,3,[],2);
    x = x(round(linspace(1,size(x,1),num_cells)),:);
    [~,order] = sort(prctile(x(:,1:100),95,2),'ascend');
    imagesc([-1/6 t_max], [1 size(x,1)],x(order,1:(t_max*12+1)),'Parent',ha(i));
    set(ha(i),'CLim',graph_lim{i},'XTick',0:4:16,'YTick',[],'LineWidth',1.5,'Box','off',...
        'FontSize',16,'TickDir','out','XTickLabel',{})
end
colormap(colormaps.byr);



% 3 x 1 bar graphs
figs.proportions(1) = figure('name','ColormapStack','Position', positionfig(520,60),'PaperPositionMode','auto');
ha = tight_subplot(1,3,[0.02 0.02]);
for i = 1:length(ha)
    freq = nfkbKO(i).fourier.freq*3600;
    off_cells = nfkbKO(i).metrics.peakfreq==0;
    osc_cells = nfkbKO(i).metrics.peakfreq>osc_cutoff;
    other_cells = (~off_cells) & (~osc_cells);
    totals = [sum(off_cells), sum(osc_cells), sum(other_cells)];
    hold(ha(i),'on')
    bar(ha(i), 3.2, sum(osc_cells)/sum(totals),0.5,'FaceColor',colors.irf,'EdgeColor','none')
    bar(ha(i), 2, sum(other_cells)/sum(totals),0.5,'FaceColor',colors.blue,'EdgeColor','none')
    hold(ha(i),'off')
    set(ha(i),'XLim',[1.2 4.8], 'YLim',[0 1], 'LineWidth',1.5,'Box','off','YTick',0:.25:1,...
        'YTickLabel',{},'XTickLabel',{},'Layer','bottom','YGrid','on','GridLineStyle',':')
end


% 3 x 1 "peakfreq"
figs2.harmonics(1) = figure('name','ColormapStack','Position', positionfig(520,100),'PaperPositionMode','auto');
ha = tight_subplot(1,3,[0.02 0.02]);
for i = 1:length(ha)
    hold(ha(i),'on')
    area(ha(i), [0.05 0.05 osc_cutoff, osc_cutoff] ,[0 0.23 0.23 0],'Facecolor',[0.741 0.918 1],'EdgeColor','none')
    area(ha(i), [osc_cutoff, osc_cutoff, 1.4 1.4] ,[0 0.23 0.23 0],'Facecolor',[0.72 1 0.72],'EdgeColor','none')
    histogram(ha(i) ,nfkbKO(i).metrics.peakfreq,nfkbKO(i).fourier.freq*3600,'Normalization','probability',...
        'FaceColor',colors.grays{3})
    hold(ha(i),'off')
    set(ha(i),'XLim',[-0.1 1.4], 'YLim',[0 0.18], 'XTick',0:0.5:2,'YTick',0:0.08:0.5,'LineWidth',1.5,'Box','on',...
        'FontSize',16,'YTickLabel',{},'XTickLabel',{},'YGrid','on','XGrid','on','Layer','top','GridLineStyle',':')
end



