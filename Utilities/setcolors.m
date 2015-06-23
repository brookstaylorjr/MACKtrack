function [colors] = setcolors(graph_flag)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% COLORS defines a selected palatte for use in graphing
%
% graph_flag    (optional) if true, show colors on heatmap
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
colors.blue = [62 107 133]/255;
colors.green = [0 195 169]/255;
colors.light_blue = [29 155 184]/255;
colors.red = [231 76 60]/255;

colors.grays = {[40 40 42]/255,...
    [74 76 79]/255,...
   [113 115 118]/255,...
   [160 165 170]/255}; 

% TLR4 model colors (Science Signaling 2015)
colors.trif = [46 122 145]/255; % blue
colors.both = [120 46 103]/255; % purple
colors.myd = [186 80 73]/255; % red
colors.irf = [48 145 50]/255;

% NFkB stimulus specificity colors
colors.tnf = [118 180 203]/255;
colors.lps = [223 79 66]/255;
colors.cpg = [136 180 69]/255;
colors.pam = [230 130 57]/255;
colors.off = [217 210 200]/255;
colors.pic = [101 76 125]/255;

if nargin>0
    if graph_flag
        names = fieldnames(colors);
        all_colors = [];
        tick_names = {};
        for i  = 1:length(names)
            if iscell(colors.(names{i}))
                for j=1:length(colors.(names{i}))
                    all_colors = cat(1,all_colors,colors.(names{i}){j});
                    tick_names = cat(1,tick_names,[names{i},'(',num2str(j),')']);
                end
            else
                all_colors = cat(1,all_colors,colors.(names{i}));
                tick_names = cat(1,tick_names,names{i});
            end
        end
        % Split colors up into groups of 4
        figure('Position',[224         913        800         600]);
        rows = ceil(length(tick_names)/4);

        h = tight_subplot(rows,1,0.1);
        disp_vect = 1:length(tick_names);

        for i = 1:rows
            axes(h(i))
            imagesc(disp_vect((i-1)*4+1:(i-1)*4+4),[1 length(tick_names)])
            set(h(i),'YTick',[])
            set(h(i),'XTick',1:4,'XTickLabel',tick_names((i-1)*4+1:(i-1)*4+4))
            colormap(all_colors)            
        end
    end
end

if nargout<1
    assignin('base','colors',colors)
end