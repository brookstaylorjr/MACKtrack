function [] = multistack_to_movie(save_name, stack_data, h_ax, time_vector)


if nargin<4
    time_vector = [];
    if nargin<3
        h_ax = gca;
    end
else
    if size(stack_data,3)~= length(time_vector)
        error('Specified time vector is incorrect length - aborting')
    end
    digits = ceil(log(max(time_vector))/log(10));
end

%%


h_fig = get(h_ax,'Parent');
cmap = get(h_fig,'Colormap');
clim = get(h_ax,'CLim');


% tmp = get(h_ax,'Children');
% full_size = size(tmp.CData);
xlim = floor(get(h_ax,'XLim')); xlim(xlim<1) = 1;
ylim = floor(get(h_ax,'YLim')); ylim(ylim<1) = 1;

z = stack_data(xlim(1):xlim(2),ylim(1):ylim(2),:);
z = (z-min(clim))/range(clim);
z(z>1) = 1; z(z<0) = 0;
z = uint8(z*255);

vidObj = VideoWriter(save_name,'Uncompressed AVI');
vidObj.FrameRate = 5;
open(vidObj)

for i = 1:size(z,3)
    R = reshape(cmap(squeeze(z(:,:,i))+1,1),[size(z,1),size(z,2)]);
    G = reshape(cmap(squeeze(z(:,:,i))+1,2),[size(z,1),size(z,2)]);
    B = reshape(cmap(squeeze(z(:,:,i))+1,3),[size(z,1),size(z,2)]);
    frame_RGB = cat(3,R,G,B);
    
    % Insert timestamp
    if ~isempty(time_vector)
        frame_RGB = insertText(frame_RGB,[8 8],[numseq(time_vector(i),digits),' min'],'TextColor','w','BoxOpacity',0,...
            'FontSize',round(range(ylim)/12),'AnchorPoint','lefttop');
    end
    % Write frame to movie
    writeVideo(vidObj,frame_RGB)
end
close(vidObj)
disp('...done.')