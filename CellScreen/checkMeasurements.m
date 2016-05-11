function [] = checkMeasurements(struct1,figtitle)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  [] = checkMeasurements(measurement_struct)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Visualize all measurements - title appropriately, show number of individuals per well, and histogram (relative)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
num_images = length(unique(struct1.image_id));
clr = linspecer(100);
clr = clr(round(linspace(1,100,num_images)),:);

num_measurements = size(struct1.measurements,2);
if nargin > 1
    figtitle = ['Diagnostic info for ',figtitle];
else
    figtitle = 'Diagnostic info for structure';
end

figure('Name',figtitle,'Position',positionfig(150*num_images,100*num_measurements)) 
ha = tight_subplot(num_measurements,num_images,[0.05 0.05], 0.1);

for i = 1:num_measurements
    bin_range = prctile(struct1.measurements(:,i),[3 97]);
    bins = linspace(bin_range(1),bin_range(2),60);
    for j= 1:num_images
        g_idx = j+(i-1)*num_images;
        histogram(struct1.measurements(struct1.image_id==j,i),bins,'FaceColor',clr(j,:),'Parent',ha(g_idx),...
            'EdgeColor','none')
        set(ha(g_idx),'XTIckLabel',{},'YTickLabel',{},'XGrid','on','XLim',bin_range,'YGrid','on')
        if j == 1
            ylabel(ha(g_idx),['Measure ',num2str(i)],'FontSize',11,'FontWeight','normal','interpreter','none')
        end
        if i ==1
            title(ha(g_idx),[struct1.image{j}(1:end-40)],'FontSize',10,'FontWeight','normal','interpreter','none')
        end
    end
end

% Sweep graphs and apply current y range
for i = 1:num_measurements
    y_max = -inf;
    y_min = inf;
    for j = 1:num_images
        g_idx = j+(i-1)*num_images; 
        if min(get(ha(g_idx),'YLim'))<y_min; y_min=min(get(ha(g_idx),'YLim'));end
        if max(get(ha(g_idx),'YLim'))>y_max; y_max=max(get(ha(g_idx),'YLim'));end
    end
    for j = 1:num_images
         g_idx = j+(i-1)*num_images;
         set(ha(g_idx),'YLim',[y_min y_max])
    end
        
        
end