function axes = diffplot(models, clusters)


%% Show 1-D progression of PPARg level
clr = cbrewer('div','RdYlBu',5);
clr = clr(end:-1:1,:);

figure('Position',positionfig(1100, 420,[3 10]));
ha = tight_subplot(3,ceil(length(all_names)/9),0.05,[0.08 0.05]);

for i = 1:size(all_cluster,2);
    hold (ha(i),'on')
    for k = 1:size(all_cluster,1)
        
        % Plot low population 1st, then high (if it exists)
        lo_idx = all_models{k,i}(:,1) == min(all_models{k,i}(:,1));
        
        on_pop = all_cluster{k,i} ~= all_lo_idx(k,i);       
        x_on = nanmedian(pparg{k}{i+starts(idx)}(on_pop));
        x_off = nanmedian(pparg{k}{i+starts(idx)}(~on_pop));
        
        
        plot(ha(i),[k-1.3 k-0.7], x_off*[1 1],'Color',colors.grays{2},...
                'LineWidth',all_models{k,i}(lo_idx,3)*18)  
                %'LineWidth',all_models{k,i}(lo_idx,3)*length(all_cluster{k,i})/1200)
        if length(lo_idx)>1
            plot(ha(i),[k-1.25 k-0.75], x_on*[1 1],'Color',colors.red,...
                'LineWidth',all_models{k,i}(~lo_idx,3)*18)
            
        end
        
        
    end
    hold (ha(i),'off')
    set(ha(i), 'YLim',[0 2500],'Box','on','Ygrid','on','XGrid','on','YTick',[0 2500 5000])
    text(mean(get(ha(i),'XLim')),max(get(ha(i),'Ylim')),[all_names{i}(5:end)],'Parent',ha(i),'Interpreter','none',...
        'HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold')
    xlabel(ha(i),'Day')
    ylabel(ha(i),'PPAR\gamma expr')

end
    legend({'Undifferentiated','Differentiated'})
