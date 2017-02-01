
function [] = checkMeasurements(struct1,figtitle)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%  [] = checkMeasurements(measurement_struct)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Visualize all measurements - title appropriately, show number of individuals per well, and histogram (relative)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
num_images = length(unique(struct1.Images));


% Get list of all measurements (of maximum size)
m_lengths = [];
m_names = fieldnames(struct1.Measurements);
for i = 1:length(m_names)
    m_lengths = cat(1,m_lengths, size(struct1.Measurements.(m_names{i}),1));
end
use_m = find( m_lengths==max(m_lengths));

if nargin > 1
    figtitle = ['Diagnostic info for ',figtitle];
else
    figtitle = 'Diagnostic info for structure';
end

clr = linspecer(100);
clr = clr(round(linspace(1,100,length(use_m))),:);


figure('Name',figtitle,'Position',positionfig(150*num_images,100*length(use_m))) 
ha = tight_subplot(length(use_m),num_images,[0.005 0.005]);

for i = 1:length(use_m)
    bin_range = prctile(struct1.Measurements.(m_names{use_m(i)}),[0.1 99.1]);
    if sum(isnan(bin_range))==0
        bins = linspace(bin_range(1),bin_range(2),60);
    else
        bins = [0 0];
        bin_range = [0 1];
        %struct1.Measurements.(m_names{use_m(i)}) = [];
    end
    for j= 1:num_images
        g_idx = j+(i-1)*num_images;
        cell_subset = struct1.CellData(:,1)==j;
        histogram(struct1.Measurements.(m_names{use_m(i)})(cell_subset),bins,'FaceColor',clr(i,:),'Parent',ha(g_idx),...
            'EdgeColor','none')
        set(ha(g_idx),'XTIckLabel',{},'YTickLabel',{},'XGrid','on','XLim',bin_range,'YGrid','on')
        if j == 1
            ylabel(ha(g_idx),m_names{use_m(i)},'FontSize',11,'FontWeight','normal','interpreter','none')
        end

        if i ==1
            title_name = struct1.Images{j};
            idx1 = strfind(title_name,'_w');
            title_name = title_name(idx1(1)-6:idx1(1)-1);
            title(ha(g_idx),title_name,'FontSize',10,'FontWeight','normal','interpreter','none')
        end
    end
end

% Sweep graphs and apply current y range
for i = 1:length(use_m)
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