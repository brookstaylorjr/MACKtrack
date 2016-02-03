function pos = positionfig(x_width, y_height,tiles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% pos = positionfig(size)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% POSITIONFIGS provides alternate figure positioning, allowing user to specify a width and
% height, then tiling figure appropriately.
% 
% INPUTS:
% x_width      desired figure width (in pixels)
% y_height     desired figure height (in pixels)
%
% OUTPUTS:
% pos          [4 x 1] figure position vector
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if nargin<3
    tiles = [2 3];
end

% Display parameters
margin = 20; % allowable distance to nearest edge of screen


% Get primary monitor position - we'll put everything on this monitor, working left-to-right
all_monitors = get(0, 'MonitorPositions');
idx = 1;
monitor_pos = all_monitors(idx,:);
while min(monitor_pos(3:4)) < 512
    try
        idx = idx+1;
        monitor_pos = all_monitors(idx,:);
    catch me
        error('All your screens are too small for figure display(< 512x512px)')
    end
end

% Divide the screen into 6 regions (either 2x3 or 3x2, depending on monitor orientation) (Leave clearance
% top and bottom, left and right)
if monitor_pos(3)>=monitor_pos(4)
    x_origins = round(linspace(margin,monitor_pos(3)-margin,tiles(2)*2 +1));
    x_origins = x_origins(2:2:end);
    y_origins = round(linspace(margin,monitor_pos(4)-margin,tiles(1)*2 +1));
    y_origins = y_origins(2:2:end);
else
    x_origins = round(linspace(margin,monitor_pos(3)-margin,tiles(1)*2 +1));
    x_origins = x_origins(2:2:end);
    y_origins = round(linspace(margin,monitor_pos(4)-margin,tiles(2)*2 +1));
    y_origins = y_origins(2:2:end);
end

[tile_y, tile_x] = meshgrid(y_origins(end:-1:1), x_origins);
tile_y = tile_y(:);tile_x = tile_x(:);


% Find most recent figure, assign it to a tile (if any) - augment one tile past that
all_figs = get(0,'Children');
if ~isempty(all_figs)
    idx = 1; val=1;
    for i = 1:length(all_figs)
        if all_figs(i).Number > val
            idx = i; val = all_figs(i).Number;
        end
    end
    
    recent_fig = all_figs(idx);
    recent_pos = get(recent_fig,'Position');
    if (recent_pos(1) > monitor_pos(3)) || (recent_pos(2) > monitor_pos(4))
        tile_idx = 1; 
    else
        recent_center = [recent_pos(1) + recent_pos(3)/2 , recent_pos(2) + recent_pos(4)/2];
        [~,~,tile_idx] = matchclosest([tile_x,tile_y],recent_center);
    end
    tile_idx  = mod(tile_idx,tiles(1)*tiles(2))+1;
else
    tile_idx  = 1;
end


% Try to center current figure in tile, but make sure we cap to upper/lower screen edges (w/ given margin)
if x_width > (monitor_pos(3)-margin*2); x_width = (monitor_pos(3)-margin*2); end
if y_height > (monitor_pos(4)-margin*2); y_height = (monitor_pos(4)-margin*2); end

fig_origin = round([tile_x(tile_idx) - x_width/2, tile_y(tile_idx) - y_height/2]);


if fig_origin(1)  < margin; fig_origin(1) = margin; end
if fig_origin(2)  < margin; fig_origin(2) = margin; end

if (fig_origin(1)+x_width) > (monitor_pos(3)-margin)
    fig_origin(1) = monitor_pos(3)-margin-x_width;
end

if (fig_origin(2)+y_height) > (monitor_pos(4)-margin)
    fig_origin(2) = monitor_pos(4)-margin-y_height;
end

pos = [fig_origin,x_width,y_height];

