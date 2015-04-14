function linePlot(measure1, CellData, options, GroupingVector, fig_handle)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Plot cell trajectories on x/y axes
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

% Make figure (if not specified)
if nargin < 5
    fig_handle = figure(gcf);% Create new figure, set properties
    set(fig_handle,'Color',[1 1 1],'InvertHardCopy', 'off','PaperPositionMode','auto')
end


% Set colormap
handles.LineColors = {[150,150,150]/255,...% gray
   [33 43 64]/255,...% dark blue
   [84,123,151]/255,...% gray-blue
   [11 89 58]/255,...% dark green 
   [186,220,221]/255}; % light blue-gray  
handles.figure1 = fig_handle;
handles.axes1 = axes('Parent', fig_handle);
handles.CellData = CellData;
handles.Measurement1 = measure1;
handles.Options = options;
handles.GroupingVector = GroupingVector;


% - - - - - - - 2) Cycle through filtered cells and plot from FrameIn to FrameOut  - - - - - - - - - - - - - - - - - - - - - - - -
handles.LineHandles = cell(size(measure1,1),1); % Cell array to hold lineseries objects
colorCycle = [1;diff(GroupingVector)];
colorInd = 1;
hold on
for i = 1:size(measure1,1)
    if (colorCycle(i)==1)
        colorInd=colorInd+1;
    end
    handles.LineHandles{i} = plot(handles.axes1, options.Times,measure1(i,:));
    set(handles.LineHandles{i},'Color',handles.LineColors{mod(colorInd,5)+1},'LineWidth',2);
end
hold off

% - - - - - - - 3) Set axis, ticks, labels, and data cursor callback- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

axis(handles.axes1,[options.TimeBounds, options.MeasurementBounds]);
set(handles.axes1,'XTick',options.TimeTicks,'YTick',options.MeasurementTicks)
set(handles.axes1,'YTickLabel',options.MeasurementTickLabels,'TickLength',[0.005 0.005])

xlabel(handles.axes1,'Time (h)','FontSize',14);
ylabel(handles.axes1,[options.Name],'FontSize',14);

dcm_obj = datacursormode(handles.figure1);
set(handles.figure1,'ResizeFcn',{@fig_resize,handles})
set(dcm_obj,'UpdateFcn',{@tooltipfcn,handles},'DisplayStyle', 'window')
% ========================================================================================


function txt = tooltipfcn(~,event_obj,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Highlight selected line and list corresponding cell
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

% Get data cursor position, format time
pos = get(event_obj,'Position');
timeHr = fix(pos(1));
timeMin = round(60*(pos(1)-timeHr));
timestr = [numseq(abs(timeHr),2),':',numseq(abs(timeMin),2)];
if (timeMin<0) || (timeHr<0)
    timestr = ['-',timestr];
end

% Make selected line wider (resetting others) and match target with a cell
target = double(get(event_obj,'Target'));


colorCycle = [1;diff(handles.GroupingVector)];
colorInd = 1;
for i = 1:size(handles.LineHandles,1)
     if target==double(handles.LineHandles{i})
    	cellIndex=i;
     end
     if (colorCycle(i)==1)
        colorInd=colorInd+1;
     end
     set(handles.LineHandles{i},'Color',handles.LineColors{mod(colorInd,5)+1},'LineWidth',2);
end
set(event_obj.Target,'LineWidth',5,'Color',[1 1 0]) 

% Set text of data cursor
if handles.Options.LogCompress
    data_pt = 10.^(pos(2));
else
    data_pt = pos(2);
end
% (Columns in CellData) 
% 1) xyPosition
% 2) index in xy pos
txt = {['XY ',num2str(handles.CellData(cellIndex,1)), ' - cell ',num2str(handles.CellData(cellIndex,2))],...
    ['Time: ',timestr],...
   [handles.Options.Name,': ',num2str(round(data_pt*10)/10)]};


% ========================================================================================

function fig_resize(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Resize plot smoothly, saving room for cursor box at bottom
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
figPos = get(handles.figure1,'Position');
set(handles.axes1,'OuterPosition',[0,60/figPos(4),1,(figPos(4)-60)/figPos(4)]);
set(handles.axes1,'LooseInset',get(handles.axes1,'TightInset')+[0 0 0.02 0.02])
% ========================================================================================
