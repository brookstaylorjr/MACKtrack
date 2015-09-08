function h = colormapStack(measure1, CellData, options, fig_handle)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Make stacked-colormap plot of cells
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Make figure (if not specified)
if nargin < 4
    fig_handle = figure(gcf);% Create new figure, set properties
    set(fig_handle,'Color',[1 1 1],'InvertHardCopy', 'off','PaperPositionMode','auto')
end


% Set figure/axes handles, attach measure1/CellData/Options to handles
handles.figure1 = fig_handle;
handles.axes1 = axes('Parent', fig_handle);
handles.CellData = CellData;
handles.Measurement1 = measure1;
handles.Options = options;


% Modify colormap so that time values outside of cell lifetime display as gray
mod_colormap = divergingmap(0:1/1023:1,[14 28 77]/255,[158 27 2]/255);
mod_colormap(1,:) = [0.1 0.1 0.1];

% Make plot, setting axes and other options
measure1DOWN = [zeros(1,size(measure1,2));measure1(1:end-1,:)];
measure1UP = [measure1(2:end,:);zeros(1,size(measure1,2))];
handles.h0 = imagesc(options.TimeBounds, [1 size(measure1,1)],nan(size(measure1)),options.MeasurementBounds);
hold on
handles.h1 = imagesc(options.TimeBounds, [1 size(measure1,1)],measure1UP,options.MeasurementBounds);
handles.h2 = imagesc(options.TimeBounds, [1 size(measure1,1)],measure1DOWN,options.MeasurementBounds);
handles.h3 = imagesc(options.TimeBounds, [1 size(measure1,1)],measure1,options.MeasurementBounds); 
hold off
colormap(mod_colormap)
c = colorbar('YTick',options.MeasurementTicks,'YTickLabel',options.MeasurementTickLabels,...
    'TickLength',[0.003 0.003]);
ylabel(c,[options.Name],'FontSize',14);set(handles.h3,'Parent',handles.axes1)
set(handles.axes1,'YTick',[],'XTick',options.TimeTicks,'TickLength',[0.005 0.005])
% - - - - - - - 4) Set fonts, labels, and data cursor callback - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
xlabel(handles.axes1,'Time (h)','FontSize',14);

dcm_obj = datacursormode(fig_handle);
set(dcm_obj,'UpdateFcn',{@tooltipfcn,handles},'DisplayStyle', 'window')
set(fig_handle,'ResizeFcn',{@fig_resize,handles},'DockControls','on')
if nargout>0
    h = handles.figure1;
end
% ========================================================================================


function txt = tooltipfcn(~,event_obj,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Highlight selected line and list corresponding cell
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

% Get data cursor position, convert and format time
pos = get(event_obj,'Position');
timeHr = fix(pos(1));
timeMin = round(60*(pos(1)-timeHr));
timestr = [numseq(abs(timeHr),2),':',numseq(abs(timeMin),2)];
if (timeMin<0) || (timeHr<0)
    timestr = ['-',timestr];
end
% Convert cursor position to image position
row = pos(2);

col_zeroed = pos(1) - handles.Options.TimeBounds(1);
col = round((col_zeroed*(size(handles.Measurement1,2)-1)/ handles.Options.TimeBounds(2))+1);


% "Magnification": use transparency to highlight selected line.
alphaMatrix1 = ones(size(handles.Measurement1));
alphaMatrix2 = ones(size(handles.Measurement1));
alphaMatrix3 = ones(size(handles.Measurement1));

backgroundRows = [row-4, row-3,row-2,row+2,row+3, row+4];
backgroundRows((backgroundRows<1) | (backgroundRows>size(handles.Measurement1,1))) = [];

upRow = row-1;
upRow(upRow<1) = [];

downRow = row+1;
downRow(downRow>size(handles.Measurement1,1)) = [];

alphaMatrix1(backgroundRows,:) = 0;
alphaMatrix2([backgroundRows,upRow],:) = 0;
alphaMatrix3([backgroundRows,upRow,downRow],:) = 0.1;

set(handles.h1,'AlphaData',alphaMatrix1)
set(handles.h2,'AlphaData',alphaMatrix2)
set(handles.h3,'AlphaData',alphaMatrix3)

% Set cursor text
if handles.Options.LogCompress
    data_pt = 10.^(handles.Measurement1(row,col));
else
    data_pt = handles.Measurement1(row,col);
end

% (Columns in CellData) 
% 1) xyPosition
% 2) index in xy pos
if ~isnan(handles.Measurement1(row,col))   
    txt = {['XY ',num2str(handles.CellData(row,1)), ' - cell ',num2str(handles.CellData(row,2))],...
    ['Time: ',timestr],...
   [handles.Options.Name,': ',num2str(round(100*data_pt)/100)]};
else
    txt = ['XY ',num2str(handles.CellData(pos(2),1)), ' - cell ',num2str(handles.CellData(pos(2),2))];

end
% ========================================================================================


function fig_resize(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Resize colorbar and figure smoothly, saving room for cursor box at bottom
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
figPos = get(handles.figure1,'Position');
set(handles.axes1,'OuterPosition',[10/figPos(3),60/figPos(4),1-50/figPos(3),(figPos(4)-60)/figPos(4)]);
% Leave room for the colorbar, but eliminate excessive margins
set(handles.axes1,'LooseInset',get(handles.axes1,'TightInset')+ [0 0 100/figPos(4) 15/figPos(4)])

% ========================================================================================
