%  TLR4 WORK , Spring 2014. Focus is on NFkB response filters (some small  component is 
% Sets to use: base dose response
% #16 - 500pg/mL LPS
% #53 - 1ng/mL LPS
% #54 - 5ng/mL LPS
% #55 - 20ng/mL LPS
% #56 - 35ng/mL LPS
% #15 - 500ng/mL LPS
% #69 - 1ug/mL LPS (+0.4ng/mL IL-4, non-physiological/no effect)
% #14 - 5ug/mL LPS

% Set environment
clear java
clc; close all
if ismac; basedir = '/Users/Brooks/Dropbox/'; %OSX
else basedir = '/home/brooks/Dropbox/'; end
savedir = [basedir,'Writing/TLR4/figs/'];
    
% Load NFkB trajectories
nfkbname = [basedir,'Code/UCSDcellTrack/Visualization/nfkb.mat'];
if ~exist(nfkbname,'file')
    ids_5ug = [16, 53, 54, 55, 56, 15, 69, 14];
    doses  = [.5 1 5 20 50 500 1000 5000];
    nfkb = struct;
    for i =1:length(ids_5ug)
        [graph] = see_nfkb(ids_5ug(i),0);
        nfkb(i).data = nan(size(graph.var));
        %  % Omit last 45min of frames
        for rw = 1:size(graph.var,1)
            inFrames = graph.celldata(rw,3):min((graph.celldata(rw,4)-9),size(graph.var,2));
            nfkb(i).data(rw,inFrames) = graph.var(rw,inFrames);
        end
        nfkb(i).data = graph.var;
        nfkb(i).dose = doses(i);
        nfkb(i).doses = doses;
    end
    save(nfkbname,'nfkb')
    clear graph info i ids
else
    load(nfkbname)
end
clear nfkbname
colors.green = [0 195 169]/255;
colors.blue = [29 155 184]/255;
colors.dark_blue = [62 107 133]/255;


%% Integral scaling in single cells
integrals = cell(size(nfkb));
for i = 1:numel(integrals)
    metrics = nfkbMetrics(nfkb(i).data,60,1);
    integrals{i} = metrics.integral;
end                

doses = log10(nfkb(1).doses);
width = 0.5;
show_bins = 0;
y_range = [-0.15 80];
bin_scale = 0.8;

% Make violin-based graph from a distribution. Bins are defined beforehand, but can be adjusted.    
vio_colors = cell(size(integrals));
for i = 1:length(integrals)
    vio_colors{i} = colors.dark_blue;
end
spaceviolin(integrals,doses,vio_colors,y_range,show_bins,width,bin_scale)
set(gca,'YTickLabel',{[]},'XTickLabel',{[]})
%print(gcf,[savedir,'subfig4_integral.eps'], '-depsc')

%%

% p65: do 1st 2 hrs vs 2-6 hrs. 
graph = see_nfkb(34,0);
metrics = nfkbMetrics(graph.var,72,0);
data_range = 0:.8:16;
metrics.vab(abs(metrics.vab)<0.1) = 0;
integral1 = nansum(metrics.vab(:,1:20),2)/24*12;
integral2 =  nansum(metrics.vab(:,21:72),2)/48*12;

% Model Gaussian distrituions
dists{1} = gmdistribution.fit(integral1,1);
dists{2}= gmdistribution.fit(integral2,2);
pdf_colors = {colors.blue/1.6, colors.green/1.8};
% Get relative frequencies
n = [hist(integral1,data_range)', hist(integral2,data_range)'];
for i = 1:size(n,2)
    n(:,i) = n(:,i)./sum(n(:,i));
end
% Draw histogram
colors = {colors.blue,colors.green};
group_width = 0.9 ; % Almost touching, 1 is touching
n_bars = size(n,2);
bar_width = group_width/n_bars;
x = data_range;
figure
hold on
shift = (group_width/(n_bars*2):group_width/n_bars:1-group_width/(n_bars*2))-.5;
for i =1:size(n,2)
    bar(x+shift(i),n(:,i),'FaceColor',colors{i},'BarWidth',bar_width,'LineStyle','none') ;
end
% Overlay modeled Gaussian Populations
pdf_shifts = [-.5 0.5];
for j = 1:length(dists)
    shift2 = pdf_shifts(j);
    dist_graph = dists{j};
    pdf_x = (min(data_range) : (max(data_range)-min(data_range))/999: max(data_range))';
    pdf_y = pdf(dist_graph,pdf_x);
    stepmult = (data_range(2)-data_range(1))/(pdf_x(2)-pdf_x(1));
    pdf_y = pdf_y/sum(pdf_y)*stepmult;
    % Solid line: overall pdf
    plot(shift2+pdf_x,pdf_y,'Color',pdf_colors{j},'LineWidth',2);
    if dist_graph.NComponents>1
        % Dashed lines: PDFs of each subpopulation
        pdf_subpop = cell(1,2);
        for i = 1:dist_graph.NComponents
            pdf_subpop{i} = normpdf(pdf_x,dist_graph.mu(i),sqrt(dist_graph.Sigma(:,:,i)));
            pdf_subpop{i} = pdf_subpop{i}/sum(pdf_subpop{i})*stepmult*dist_graph.PComponents(i);
            plot(pdf_x+shift2,pdf_subpop{i},'--','Color',pdf_colors{j},'LineWidth',2);
        end   
    end

end
hold off


%% - - - - - - - - - METHODOLOGY - - - - - - - - - 
% a) Normalizing NFkB signal  
id = 15;

[graph,info,measure] = see_nfkb(id,0);
nfkb_cyto = measure.NFkBCytoplasm(info.keep,:);
nfkb_cyto = nfkb_cyto./repmat(nanmean(nfkb_cyto(:,1:3),2),1,size(nfkb_cyto,2));

nfkb_nuc = measure.NFkBNuclear(info.keep,:);
nfkb_nuc = nfkb_nuc./repmat(nanmean(nfkb_nuc(:,1:3),2),1,size(nfkb_nuc,2));


nfkb_corr = graph.var;

t = 0:1/12:(size(nfkb_nuc,2)-1)/12;
t2 = 0:1/12:(size(nfkb_corr,2)-1)/12;


mod_colormap = divergingmap(0:1/1023:1,[12 12 77]/255,[158 4 0]/255);

figure(1),imagesc(nfkb_cyto), colormap(mod_colormap)
figure(2), imagesc(nfkb_nuc), colormap(mod_colormap)

figure(3)
plot(t, nanmean(nfkb_nuc),'LineWidth',2,'Color',colors.blue)
hold on
plot(t,nanmean(nfkb_cyto),'LineWidth',2,'Color',colors.green)
plot(t2,nanmean(nfkb_corr),'LineWidth',2,'Color',colors.red)
hold off
%%

ordr = randperm(size(nfkb_corr,1),5);
figure(4)
plot(t2',nfkb_corr(ordr,:)','-','LineWidth',2)
hold on
plot(t',nfkb_nuc(ordr,:)','.','LineWidth',2)
hold off

%%

data = measure.NFkBNuclear(info.keep,1:144);
local_std = stdfilt(data,ones(1,7));
target_windows = local_std<repmat(prctile(local_std,25,2),1,size(local_std,2));
tmp = data;
tmp(~target_windows) = nan;
tmp(:, 1:24) = nan;
base_vals = prctile(tmp,25,2);
start_vals = nanmean(data(:,1:3),2);

for i =1:size(data,1)
    data(i,1:24) = data(i,1:24) - linspace(start_vals(i)-base_vals(i),0,24);
end

data = data./repmat(nanmean(data(:,1:3),2),1,size(data,2));

ordr = randperm(size(data,1),8);
t = 0:1/12:(size(data,2)-1)/12;
plot(t,data(ordr,:)','.')
%%

figure,imagesc(stdfilt(data,ones(1,7)),[0 2000])





%% b) Behavior vs RelA brightness
id = 14;
[graph,info,measure] = see_nfkb(id,0);

metrics = nfkbMetrics(graph.var,96,1);
colors = setcolors;
% Total integrated activity (8 hrs) vs brightness
x_var = prctile(measure.NFkBNuclear(info.keep,1:8),18.75,2);
y_var = metrics.integral;
r = corrcoef(x_var,y_var);
figure,plot(x_var,y_var,'.','MarkerSize',20, 'Color',colors.green), title(['r = ',num2str(round(r(1,2)*100)/100)])
xlabel('Starting nuclear RelA (a.u.)')
ylabel('Integrated translocation activity')
axis([0 40 0 190])
% print(gcf,[basedir,'brightness_vs_activity(pk1).eps'], '-depsc')

% First-peak amplitude vs brightness
x_var = prctile(measure.NFkBNuclear(info.keep,1:8),18.75,2);
y_var = metrics.pk1_integral;
drops = isnan(x_var)|isnan(y_var);
r = corrcoef(x_var(~drops),y_var(~drops));
figure,plot(x_var,y_var,'.','MarkerSize',20, 'Color',colors.green), title(['r = ',num2str(round(r(1,2)*100)/100)])
xlabel('Starting nuclear RelA (a.u.)')
ylabel('First peak integral')
axis([0 40 0 40])
% print(gcf,[basedir,'brightness_vs_activity(pk1).eps'], '-depsc')

%% 
num_bins = 40;
x_range = [0.-.1 2.5];
metrics = nfkbMetrics(nfkb(ind).data,120,0);
for ind = 1:9
    mod_colormap = divergingmap(0:1/(num_bins-1):1,[12 12 77]/255,[158 4 0]/255);
    x = min(x_range):diff(x_range)/(num_bins-1):max(x_range);
    all_hists = zeros(num_bins,80);
    for i = 1:80
           n = hist(metrics.vab(:,i),x);
           all_hists(end:-1:1,i) = n./sum(n);
    end

    figure,imagesc(all_hists,[0 .15]),colormap(mod_colormap)
end
set(gcf,'Position',  [283   824   908   134])






% The base TLR response is oscillatory, but single cells are desynchronized

% a) Fourier signatures, by subpopulation (- 1/f noise)
id = 9;
var1 = nfkb(id).data;
num_clst = 1;

graph.t = 0:1/12:144/12;

metrics = nfkbMetrics(var1,30,1);
Gaussian.Dist = gmdistribution.fit(metrics.pk1_integral,num_clst); % 2 subpopulations
% Order Gaussian mixture components by mean and cluster
[~,sortind] = sort(Gaussian.Dist.mu,'ascend');
[~, sortorder] = sort(sortind,'ascend');
Gaussian.Cluster = cluster(Gaussian.Dist,metrics.integral);
tmpcol = Gaussian.Cluster;
tmpcol(~isnan(tmpcol)) = sortorder(tmpcol(~isnan(tmpcol)));
Gaussian.Cluster = tmpcol;

colors1 = [colors.blue; colors.green];
for clst = 1:num_clst
    subpop1 = metrics.vab(Gaussian.Cluster==clst,:);
    idx = 2*(clst-1)+1;
    figure
    set(gcf,'Position',[563   907   984   427],'PaperPositionMode','auto')
    hold on
    for i = 1:size(subpop1,1)
        fill([graph.t,graph.t(end:-1:1)],[subpop1(i,:),zeros(1,size(subpop1,2))],colors1(clst,:)/(1+i/100),...
            'FaceAlpha',0.025,'EdgeColor','none')
    end
    hold off
    set(gca,'XLim',[0 12],'YLim',[-.01 3])
    set(gca,'Visible','off')
    %print(gcf,[basedir,nfkb(id).name,'-traces(clst',num2str(clst),').eps'], '-depsc')

    % Bar graph
    figure
    set(gcf,'Position',[ 418         904        1094         431],'PaperPositionMode','auto')
    % Take a guess at 1/f noise based on 2nd value...

    % NOTE: may be better to do this on a cell-by-cell basis. Try both!!!

    avg_fft = mean(metrics.fourier(Gaussian.Cluster==clst,:));
    flicker = 1./metrics.freq*avg_fft(2)*metrics.freq(2);

    bar(metrics.freq,avg_fft-flicker,...
        'FaceColor',colors1(clst,:),'BarWidth',0.9,'LineStyle','none')
    set(gca,'YTickLabel',{[]})
    ylabel('|H(f)|')
    xlabel('Frequency (Hz)')
    set(gca,'XLim',[0 8e-4],'YLim',[0 0.2])
   %print(gcf,[basedir,nfkb(id).name,'-fourier(clst',num2str(clst),').eps'], '-depsc')
end

%% Oscillations: behavior vs brightness
id = 70;
[graph,info,measure] = see_nfkb(id,0);

metrics = nfkbMetrics(graph.var,100,1);
brightness = nanmedian(measure.NFkBCytoplasm(info.keep,1:3),2);
areas = nanmean(measure.Area(info.keep,1:3),2);

x_var = brightness;
y_var = metrics.pk1_integral;

r = corrcoef(x_var,y_var);
figure,plot(x_var,y_var,'.','MarkerSize',20, 'Color',colors.green), title(['r = ',num2str(round(r(1,2)*100)/100)])
xlabel('Cytoplasmic RelA (a.u.)')
ylabel('First peak translocation activity')
axis([0 45 0 40])
print(gcf,[basedir,'brightness_vs_activity(pk1).eps'], '-depsc')

x_var = brightness;
y_var = metrics.integral;

r = corrcoef(x_var,y_var);
figure,plot(x_var,y_var,'.','MarkerSize',20, 'Color',colors.green), title(['r = ',num2str(round(r(1,2)*100)/100)])
xlabel('Cytoplasmic RelA (a.u.)')
ylabel('Total translocation activity')
axis([0 45 0 100])
print(gcf,[basedir,'brightness_vs_activity(tot).eps'], '-depsc')




%% Looking at average behavior: dynamic features evident in population response
mod_colormap = divergingmap(0:1/(length(nfkb)-1):1,[12 12 77]/255,[158 4 0]/255);
figure
for i =1:length(nfkb)
    %metrics =  nfkbMetrics(nfkb(i).data,120,0);
    mean_activation = nanmean(nfkb(i).data);
    plot(0:1/12:(12*12)/12, mean_activation(1:145),'Color',mod_colormap(i,:),'LineWidth',2)
    hold on
end
hold off
xlabel('Time'),ylabel('Mean value above baseline')

figure
for i =1:length(nfkb)
    %metrics =  nfkbMetrics(nfkb(i).data,120,0);
     mean_activation = nanmean(nfkb(i).data);
    plot(0:1/12:51/12, mean_activation(3:54),'Color',mod_colormap(i,:),'LineWidth',2)
    hold on

    if (i ==7) || (i==3)
        plot([0 .5 1 2 4], mean_activation([3, 9, 15, 27, 50]),'.','MarkerSize',40,'Color',mod_colormap(i,:))
    end
end
hold off
xlabel('Time'),ylabel('Mean value above baseline')
axis([0 4.5 0 1.5])


%% - - - - Summary plot: 1st peak vs later activity - - - -
doses = nfkb(1).doses;
idxs = 1:length(nfkb);
mus = zeros(length(nfkb),2);
stds = zeros(length(nfkb),2);

for i = 1:length(nfkb)
    idx = idxs(i);
    metrics = nfkbMetrics(nfkb(idx).data,80,1);
    mus(i,1) = mean(metrics.pk1_integral);
    stds(i,1) = std(metrics.pk1_integral);
    mus(i,2) = mean(metrics.integral);
    stds(i,2) = std(metrics.integral);

end

stds(:,1) = stds(:,1);
stds(:,2) = stds(:,2);
x = 1;

figure
hold on;
set(gca,'ColorOrder',[colors.blue;colors.green]/x,'XTick',[0 2 4],'XTickLabel',{'1','100','10000'})
plot(log10(doses),mus,'.','MarkerSize',25);
errorbar(log10(doses),mus(:,1),stds(:,1),'--','Color',colors.blue/x,'LineWidth',2);
errorbar(log10(doses),mus(:,2),stds(:,2),'--','Color',colors.green/x,'LineWidth',2);
hold off
axis([-1 4 -1 60])
legend({'Primary activity','Total activity'})
%print(gcf,['/home/brooks/Dropbox/Presentations/ps+ai/2014-03_labmtg/primary-secondary.eps'], '-depsc')


%%  IFNg vs control: oscillation dropout?
clear java; clc;

if ismac; basedir = '/Users/Brooks/Dropbox/'; %OSX
else basedir = '/home/brooks/Dropbox/'; end

% Load NFkB trajectories
filename = [basedir,'Code/UCSDcellTrack/Visualization/nfkb_M1M2.mat'];
if ~exist(filename,'file')
    ids = [69 65 72];
    doses  = [1000 1000 1000];
    nfkb = struct;
    for i =1:length(ids)
        [graph] = see_nfkb(ids(i),0);
        nfkb(i).data = nan(size(graph.var));
        %  % Omit last 45min of frames
        for rw = 1:size(graph.var,1)
            inFrames = graph.celldata(rw,3):min((graph.celldata(rw,4)-9),size(graph.var,2));
            nfkb(i).data(rw,inFrames) = graph.var(rw,inFrames);
        end
        nfkb(i).data = graph.var;
        nfkb(i).dose = doses(i);
        nfkb(i).doses = doses;
    end
    save(filename,'nfkb')
    clear graph info i ids
else
    load(filename)
end
clear polarname
savedir = '/Users/Brooks/Dropbox (Biodynamics Lab)/Writing/2014 TLR Grant/Figures/raw';
colors.trif = [46 122 145]/255; % blue
colors.both = [120 46 103]/255; % purple
colors.myd = [186 80 73]/255; % red
colors.irf = [48 145 50]/255;

%%
fig1 = figure('Position',[480 600 500  500], 'PaperPositionMode','auto');
duration_sum = zeros(100,length(nfkb));
top = 1.7;
t_len =  60;
vect = 0:top/99:top;
for i = 1:length(nfkb)
    nfkb1 = nfkb(i).data(:,1:t_len)+0.15;
    nfkb1(sum(isnan(nfkb1),2)>3,:) = [];
    for j = 1:length(duration_sum)
        all_durations = zeros(1,size(nfkb1,1));
        for k = 1:size(nfkb1,1)
            q = diff([0 2*(nfkb1(k,:)>vect(j)) 0] == 2);
            v = find(q == -1) - find(q == 1);
            if ~isempty(v)
                all_durations(k) = max(v)/12;
            end
        end
        duration_sum(j,i) = mean(all_durations);
    end
end
linecolors = linspecer(length(nfkb),'sequential');
set(gcf,'DefaultAxesColorOrder',linecolors)
plot(vect,duration_sum,'LineWidth',2),axis([min(vect) max(vect) 0 t_len/12])
set(gca,'XTickLabel',{},'YTickLabel',{},'LineWidth',1.5)
%print(fig1,[savedir,'/duration-M1.eps'], '-depsc') 
%%
fig1 = figure('Position',[480 600 500  500], 'PaperPositionMode','auto');
plot(vect,diff(duration_sum,[],2), 'LineWidth', 2), axis([0 max(vect) -.1 1.5])
%print(fig1,[savedir,'/duration-diff-M1.eps'], '-depsc') 
%%
nfkb = nfkb([1 2]);
fig1 = figure('Position',[480 600 500  500], 'PaperPositionMode','auto');
duration_sum = zeros(100,length(nfkb));
top = 2.0;
t_len =  120;
vect = 0:top/99:top;
for i = 1:length(nfkb)
    nfkb1 = nfkb(i).data(:,1:t_len)+0.15;
    nfkb1(sum(isnan(nfkb1),2)>3,:) = [];
    size(nfkb1)
    for j = 1:length(duration_sum)
        duration_sum(j,i) = nansum(nfkb1(:)>vect(j))/size(nfkb1,1)/12;
    end
end
linecolors = linspecer(length(nfkb),'sequential');
set(gcf,'DefaultAxesColorOrder',linecolors)
plot(vect,duration_sum,'LineWidth',2),axis([min(vect) max(vect) 0 t_len/12])
set(gca,'XTickLabel',{},'YTickLabel',{},'LineWidth',1.5)
print(fig1,[savedir,'subfig3_duration(data).eps'], '-depsc')


fig2 = figure('Position',[480 600 300  100], 'PaperPositionMode','auto');
t_len =  120;
thresh = 0.8;
duration_examp = zeros(size(nfkb));
hold on
for i = 1:length(nfkb)
    nfkb1 = nfkb(i).data(:,1:t_len)+0.15;
    nfkb1(sum(isnan(nfkb1),2)>3,:) = [];
    duration_examp(i) = nansum(nfkb1(:)>thresh)/size(nfkb1,1)/12;
end
plot(nfkb(1).doses,duration_examp, 'LineWidth',1.5,'Color',[0 0 0])
hold on
for i = 1:length(nfkb)
    plot(nfkb(1).doses(i),duration_examp(i),'o','Color',linecolors(i,:),'LineWidth',2,'MarkerFaceColor',[1 1 1],'MarkerSize',8)
end
hold off
set(gca,'XTickLabel',{},'YTickLabel',{},'LineWidth',1.5,'XScale','log',...
    'XLim',[0.25 6000],'YLim',[0 2],'XTick',[.5 50 5000],'Box','on')

%%
% Graphing parameters
n_cells = 160; % # of cells displayed (overlay)
wobble = 80; % Color variation (overlay)
desat_rate = 1; % Speed of color desaturation
sing_cells = 6; % Number of single cells to graph

% MAIN LOOP
% Create data/ data-specific params
t = 0:1/12:6; % (hrs)

color{1} = [6 133 135]; 
color{2} = [17 47 65]; 

sets = [2 1];

y_range = [-0.05 2.7];
max_alpha = 4;
data = cell(1,numel(sets));
for i = 1:numel(sets)
    data{i} = nfkb(sets(i)).data(:,1:length(t));
end      


% - - - 1) Transparency overlays - - - 
desat = desat_rate.^(linspace(0,1,numel(data))); % Desaturation function
h1 = figure('Position',[37, 50, 300, 120*numel(data)],'PaperPositionMode','auto');
ha = tight_subplot(numel(data),1,0);
h2 = figure('Position',[37, 50, 300, 120*numel(data)],'PaperPositionMode','auto','Visible','off');
ha_lines = tight_subplot(numel(data),1,0);
ticks = maketicks(t, y_range,0);

for ind = 1:numel(data)
    axes(ha(ind)), hold on
    ind2 = numel(data)+1-ind;
    rows = randperm(size(data{ind2},1));
    rows = rows(1:n_cells);
    for i = 1:length(rows)
        tmp = wobble*(rand(1,3)-[.5 .5 .5])/255;
        fill_color = (color{ind}/255/desat(ind))+tmp;
        fill_color(fill_color<0) = 0;
        fill_color(fill_color>1) = 1;        
        fill([t,t(end:-1:1)],...
            [data{ind2}(rows(i),:),y_range(1)*ones(1,size(data{ind2},2))],...
            fill_color,'FaceAlpha',max_alpha/n_cells,'EdgeColor','none')
    end
    hold off
    set(ha(ind),'XTickLabel',{[]})
    set(ha(ind),'YTick',[],'YTickLabel',{[]}, 'XTick',[])
    set(ha(ind),'XLim',[min(t) max(t)],'YLim', y_range)
    set(ha_lines(ind),'XLim',[min(t) max(t)],'YLim', y_range)
    set(ha_lines(ind),'XTick',ticks.TimeTicks,'YTick',[], 'XTickLabel',{[]},'Box','on','LineWidth',1.5)
end

print(h1,[savedir,'/overlay-M1.tiff'], '-dtiff')
print(h2,[savedir,'/overlaylines-M1.eps'], '-deps') 
    
    
%%
t = 0:1/12:(size(nfkb(1).data,2)-1)/12;
avg_ctrl = nanmean(nfkb(1).data);
avg_m1 = nanmean(nfkb(2).data);
fig1 = figure('Position',[490 1066   533   294],'PaperPositionMode','auto');
plot(t(1:end-2),avg_m1(1:end-2),t(1:end-2),avg_ctrl(3:end),'LineWidth',2),axis([0 10 -0.05 1.8])
set(gca,'XTickLabel',{},'YTickLabel',{},'YTick', [0 0.5 1 1.5],'XTick',[0 5 10],'LineWidth',1.5)
print(fig1,[savedir,'/avg-M1.eps'], '-depsc') 

fig1 = figure('Position',[490 1066   533   150],'PaperPositionMode','auto');
plot(t(1:end-2),(avg_m1(1:end-2)-avg_ctrl(3:end))./max(avg_ctrl),'LineWidth',2,'Color',[.8 .8 .8]), axis([0 10 -.5 .5])
set(gca,'XTickLabel',{},'YTickLabel',{},'YTick', [-.5 0 .5],'XTick',[0 5 10],'LineWidth',1.5)

print(fig1,[savedir,'/avg-diff-M1.eps'], '-depsc') 


%% LOAD+SAVE BLOCKER EXPERIMENTS
clear java; clc;

if ismac; basedir = '/Users/Brooks/Dropbox/'; %OSX
else basedir = '/home/brooks/Dropbox/'; end

% Load NFkB trajectories
filename = [basedir,'Code/UCSDcellTrack/Visualization/nfkb_blockers.mat'];
if ~exist(filename,'file')
    ids = 33:36;
    doses  = [500 500 500 500];
    nfkb = struct;
    for i =1:length(ids)
        [graph] = see_nfkb(ids(i),0);
        nfkb(i).data = nan(size(graph.var));
        %  % Omit last 45min of frames
        for rw = 1:size(graph.var,1)
            inFrames = graph.celldata(rw,3):min((graph.celldata(rw,4)-9),size(graph.var,2));
            nfkb(i).data(rw,inFrames) = graph.var(rw,inFrames);
        end
        nfkb(i).data = graph.var;
        nfkb(i).dose = doses(i);
        nfkb(i).doses = doses;
    end
    save(filename,'nfkb')
    clear graph info i ids
else
    load(filename)
end
clear polarname
savedir = '/home/brooks/Dropbox (Biodynamics Lab)/Writing/2014 TLR Grant/Figures/raw';
colors.trif = [46 122 145]/255; % blue
colors.both = [120 46 103]/255; % purple
colors.myd = [186 80 73]/255; % red
colors.irf = [48 145 50]/255;



%% TOTAL Duration: ctrl dose response, plus M1-polarized
load([basedir,'Code/UCSDcellTrack/Visualization/nfkb.mat']);
findval = @(val,vect) find(abs(vect-val)==min(abs(vect-val)),1,'first');

thresh1 = 0.35;


fig1 = figure('Position',[480 600 500  500], 'PaperPositionMode','auto');
top = 1.6;
t_len =  120;
vect = 0:top/99:top;
sets = [1 4 6 7 8];
duration_sum = zeros(100,length(sets));

for i = 1:length(sets)
    nfkb1 = nfkb(sets(i)).data(:,1:t_len)+0.15;
    nfkb1(sum(isnan(nfkb1),2)>3,:) = [];
    size(nfkb1)
    for j = 1:length(duration_sum)
        duration_sum(j,i) = nansum(nfkb1(:)>vect(j))/size(nfkb1,1)/12;
    end
end
linecolors = linspecer(length(sets),'sequential');
set(gcf,'DefaultAxesColorOrder',linecolors)
plot(vect,duration_sum,'LineWidth',2),axis([min(vect) max(vect) 0 t_len/12])
set(gca,'XTickLabel',{},'YTickLabel',{},'LineWidth',1.5,'XTick',[0 0.5 1.0 1.5],'YTick',[0 2 4 6 8 10])

sums = zeros(length(sets)+1,1);
sums(1:length(sets)) = duration_sum(findval(thresh1,vect),:);

% Now load M1 set/display it
load( [basedir,'Code/UCSDcellTrack/Visualization/nfkb_M1M2.mat'])
duration_sum = zeros(100,1);
for i = 2
    nfkb1 = nfkb(i).data(:,1:t_len)+0.15;
    nfkb1(sum(isnan(nfkb1),2)>3,:) = [];
    size(nfkb1)
    for j = 1:length(duration_sum)
        duration_sum(j) = nansum(nfkb1(:)>vect(j))/size(nfkb1,1)/12;
    end
end
hold on
plot(vect,duration_sum,'LineWidth',2,'Color',[0.2 0.2 0.2])
hold off
sums(end) = duration_sum(findval(thresh1,vect),:);



figure('PaperPositionMode','auto','Position',[  500   989   434   361])
bar(1:5,sums(1:end-1),'FaceColor',[6 133 135]/255,'EdgeColor','none')
hold on
bar(7.2,sums(end),'FaceColor',[27 52 95]/255)
hold off
set(gca,'XTickLabel',{},'YTick', 0:2:6,'XLim',[0 8.2],'YTickLabel',{},'LineWidth',1.5,'XTick',[])

% print(fig1,[savedir,'/total-duration-M1.eps'], '-depsc') 


% CONSECUTIVE Duration: ctrl dose response, plus M1-polarized
findval = @(val,vect) find(abs(vect-val)==min(abs(vect-val)),1,'first');

load([basedir,'Code/UCSDcellTrack/Visualization/nfkb.mat']);
fig1 = figure('Position',[480 600 500  500], 'PaperPositionMode','auto');
top = 1.6;
t_len =  60;
vect = 0:top/99:top;
duration_sum = zeros(100,length(sets));

for i = 1:length(sets)
    nfkb1 = nfkb(sets(i)).data(:,1:t_len)+0.15;
    nfkb1(sum(isnan(nfkb1),2)>3,:) = [];
    for j = 1:length(duration_sum)
        all_durations = zeros(1,size(nfkb1,1));
        for k = 1:size(nfkb1,1)
            q = diff([0 2*(nfkb1(k,:)>vect(j)) 0] == 2);
            v = find(q == -1) - find(q == 1);
            if ~isempty(v)
                all_durations(k) = max(v)/12;
            end
        end
        duration_sum(j,i) = mean(all_durations);
    end
end
linecolors = linspecer(length(sets),'sequential');
set(gcf,'DefaultAxesColorOrder',linecolors)
plot(vect,duration_sum,'LineWidth',2),axis([min(vect) max(vect) 0 t_len/12])
set(gca,'XTickLabel',{},'YTickLabel',{},'LineWidth',1.5,'XTick',[0 0.5 1.0 1.5],'YTick',[0 1 2 3 4 5])


sums = zeros(length(sets)+1,1);
sums(1:length(sets)) = duration_sum(findval(thresh1,vect),:);

% Now load M1 set/display it
load( [basedir,'Code/UCSDcellTrack/Visualization/nfkb_M1M2.mat'])
duration_sum = zeros(100,1);

for i = 2
    nfkb1 = nfkb(i).data(:,1:t_len)+0.15;
    nfkb1(sum(isnan(nfkb1),2)>3,:) = [];
    for j = 1:length(duration_sum)
        all_durations = zeros(1,size(nfkb1,1));
        for k = 1:size(nfkb1,1)
            q = diff([0 2*(nfkb1(k,:)>vect(j)) 0] == 2);
            v = find(q == -1) - find(q == 1);
            if ~isempty(v)
                all_durations(k) = max(v)/12;
            end
        end
        duration_sum(j) = mean(all_durations);
    end
end
hold on
plot(vect,duration_sum,'LineWidth',2,'Color',[0.2 0.2 0.2])
hold off
sums(end) = duration_sum(findval(thresh1,vect),:);
% print(fig1,[savedir,'/consecutive-duration-M1.eps'], '-depsc') 

figure('PaperPositionMode','auto','Position',[  500   989   434   361])
bar(1:5,sums(1:end-1),'FaceColor',[6 133 135]/255,'EdgeColor','none')
hold on
bar(7.2,sums(end),'FaceColor',[27 52 95]/255)
hold off
set(gca,'XTickLabel',{},'YTick', 0:3,'XLim',[0 8.2],'YTickLabel',{},'LineWidth',1.5,'XTick',[])
%% Load/show colormap of  CpG data

% Load the CpG data (it's on my mac desktop)
i = 1;
ids = [69 65 72];
doses  = [10];
nfkb = struct;
[graph] = see_nfkb_old('/Users/Brooks/Desktop/AllMeasurements.mat',0);
nfkb(i).data = nan(size(graph.var));
%  % Omit last 45min of frames
for rw = 1:size(graph.var,1)
    inFrames = graph.celldata(rw,3):min((graph.celldata(rw,4)-9),size(graph.var,2));
    nfkb(i).data(rw,inFrames) = graph.var(rw,inFrames);
end
nfkb(i).data = graph.var;
nfkb(i).dose = 10;
nfkb(i).doses = doses;

%%
doses = 1;
mod_colormap = divergingmap(0:1/1023:1,[12 12 77]/255,[158 4 0]/255);
mod_colormap(1,:) = [0.1 0.1 0.1];
fig1 = figure('Position',[490   710   402   200],'PaperPositionMode','auto');
ha = tight_subplot(1,1);colormap(mod_colormap)

for i =1:length(doses)
    data1 = nfkb(doses(i)).data;
    order1 = hierarchial(data1,0);
    data1 = data1(order1,:);
    imagesc(data1(round(linspace(1,size(data1,1),200)),1:109),'Parent',ha(i))
    set(ha(i),'CLim',[-0.1 2],'XTick',[1 37 73 109 145],'YTick',[],'XTickLabel',{},'Box','off','TickDir','out') 
    if i ==length(doses)
           set(ha(i),'XTickLabel',{'0', '3', '6', '9', '12'})
    end
end
savedir = '/Users/Brooks/Dropbox (Biodynamics Lab)/Writing/2014 TLR Grant/Figures/raw';
print(fig1,[savedir,'/CpG-response.eps'], '-depsc') 
