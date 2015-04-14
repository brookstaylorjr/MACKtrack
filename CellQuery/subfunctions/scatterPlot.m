function scatterPlot(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SCATTERPLOT Display scatter plot from handles.Export data
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Initialize axes and slider and set properties
if isempty(get(handles.PlotFigure,'Children'))
    % Axes1: measurement 1
    handles.axes1 = axes('Parent',handles.PlotFigure);
    set(handles.axes1,'LooseInset',[0.13, 0.11, 0.095, 0.075]/2)
    % Initialize slider
    stepSize = 1/size(handles.Export.Measurement1,2);
    handles.slider1  = uicontrol('Style', 'slider','Max',size(handles.Export.Measurement1,2),...
        'Min',1,'Value',1,'SliderStep', [stepSize stepSize*10]);
    set(handles.slider1,'Callback',{@slider_Callback,handles},'BackgroundColor',[.9 .9 .9])
    % Initialize "save" button
    handles.pushbutton1  = uicontrol('Style', 'pushbutton','String','Save Series');
    set(handles.pushbutton1,'Callback',{@pushbutton_Callback,handles},'BackgroundColor',[.9 .9 .9])
    % Figure properties
    set(handles.PlotFigure,'ResizeFcn',{@fig_resize,handles})
end




% Set current frame, then display
frame = get(handles.slider1,'Value');
formplot(handles,frame)


% ========================================================================================


function formplot(handles, frame)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Set colors and form scatter plot
%
% handles     master structure
% frame       current time point
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Set colors
point_colors = [[194,224,120]/255;...% yellow-green
   [33 43 64]/255;...% dark blue
   [84,123,151]/255;...% gray-blue
   [11 89 58]/255;...% dark green 
   [186,220,221]/255]; % light blue-gray

if isfield(handles.Export,'Gaussian')
    tmpcol = handles.Export.Gaussian.Cluster(:,frame);
    tmpcol(isnan(tmpcol)) = 4;
    tmpcol = mod(tmpcol,4)+1;
    colors = point_colors(tmpcol,:);
else
    colors = point_colors(mod(handles.Export.GroupingVector,4)+1,:);
end

% Make scatter plot
scatter(handles.axes1, handles.Export.Measurement1(:,frame),handles.Export.Measurement2(:,frame), 50,colors,'filled')  

% Update time element and set title
% Update time element
pos = handles.Export.Measurement1Info.Times(frame);
timeHr = fix(pos);
timeMin = round(60*(pos-timeHr));
titlestr = [numseq(abs(timeHr),2),':',numseq(abs(timeMin),2)];
if (timeMin<0) || (timeHr<0)
    titlestr = ['-',titlestr];
end
title(handles.axes1,titlestr,'FontSize',16,'FontName','Arial')

% Set axis properties
xlabel(handles.axes1,[handles.Export.Measurement1Info.Name],'FontSize',14);
ylabel(handles.axes1,[handles.Export.Measurement2Info.Name],'FontSize',14);
set(handles.axes1,'XTick',handles.Export.Measurement1Info.MeasurementTicks,'XTickLabel',handles.Export.Measurement1Info.MeasurementTickLabels,...
'YTick',handles.Export.Measurement2Info.MeasurementTicks,'YTickLabel',handles.Export.Measurement2Info.MeasurementTickLabels,...
'TickLength',[0.005 0.005]);
axis(handles.axes1,[handles.Export.Measurement1Info.MeasurementBounds,handles.Export.Measurement2Info.MeasurementBounds])
set(handles.PlotFigure,'ToolBar','figure')
dcm_obj = datacursormode(handles.PlotFigure);
set(dcm_obj,'UpdateFcn',{@tooltipfcn,handles},'DisplayStyle', 'datatip')


% Display contours for 2-D Gaussian modeling
if strcmp(handles.Options.Grouping,get(handles.radiobutton8E,'String'))    
    hold(handles.axes1,'on')
    % Solid line: PDF of combined function
    pdf_x = (handles.Export.Measurement1Info.MeasurementBounds(1) :...
        (handles.Export.Measurement1Info.MeasurementBounds(2)-handles.Export.Measurement1Info.MeasurementBounds(1))/99: ...
        handles.Export.Measurement1Info.MeasurementBounds(2))';
    pdf_y = (handles.Export.Measurement2Info.MeasurementBounds(1) :...
        (handles.Export.Measurement2Info.MeasurementBounds(2)-handles.Export.Measurement2Info.MeasurementBounds(1))/99: ...
        handles.Export.Measurement2Info.MeasurementBounds(2))';    
    [mesh1 mesh2] = meshgrid(pdf_x,pdf_y);
    pdf_z = pdf(handles.Export.Gaussian.Dist{frame},[reshape(mesh1,numel(mesh1),1),reshape(mesh2,numel(mesh2),1)]);
    pdf_z = reshape(pdf_z,numel(pdf_x),numel(pdf_y));
    pdf_z = pdf_z/sum(pdf_z(:));
    contour(handles.axes1,mesh1,mesh2,pdf_z,[.0001 .0002 .0005 .001 .005],'LineWidth',2), colormap gray
    set(handles.axes1,'CLim',[0 .008]) % Set top bound at 50x avg value
    % Dashed lines: PDFs of each subpopulation
%     pdf_subpop = cell(1,handles.Options.Subpopulations);
%     for i = 1:handles.Options.Subpopulations
%         pdf_subpop{i} = normpdf(pdf_x,handles.Export.Gaussian.Dist{frame}.mu(i,axesnum),sqrt(handles.Export.Gaussian.Dist{frame}.Sigma(axesnum,axesnum,i)));
%         pdf_subpop{i} = pdf_subpop{i}/sum(pdf_subpop{i})*stepmult*handles.Export.Gaussian.Dist{frame}.PComponents(i);
%         plot(handles.(axesname),pdf_x,pdf_subpop{i},'--','Color',[0.2 0.2 0.2],'LineWidth',2);
%     end    
    hold(handles.axes1,'off')

end
% ========================================================================================


function txt = tooltipfcn(~,event_obj,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Highlight selected line and list corresponding cell
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
% Get data cursor position, format time
pos = get(event_obj,'Position');

% Set text of data cursor
if handles.Options.LogCompress1
    data_pt_x = 10.^(pos(1));
else
    data_pt_x = pos(1);
end

if handles.Options.LogCompress2
    data_pt_y = 10.^(pos(2));
else
    data_pt_y = pos(2);
end
frame = get(handles.slider1,'Value');

cellIndex = find((handles.Export.Measurement1(:,frame)==pos(1)) & (handles.Export.Measurement2(:,frame)==pos(2)), 1);
% (Columns in CellData) 
% 1) xyPosition
% 2) index in xy pos
txt = {['XY ',num2str(handles.Export.CellData(cellIndex,1)), ' - cell ',num2str(handles.Export.CellData(cellIndex,2))],...
   [handles.Export.MeasureField1,': ',num2str(round(data_pt_x*10)/10)],...
   [handles.Export.MeasureField2,': ',num2str(round(data_pt_y*10)/10)]};
% ========================================================================================


function slider_Callback(hObject,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Update time with slider movement/ redraw graph
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
set(hObject,'Value',round(get(hObject,'Value')))
scatterPlot(handles)
% ========================================================================================


function fig_resize(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Resize scatter plot and slider smoothly
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
figPos = get(handles.PlotFigure,'Position');
h = figPos(4);
w = figPos(3);
set(handles.slider1,'Position',[w/4-50 5 w/2 20]);
set(handles.pushbutton1,'Position',[3*w/4-40 4 100 24]); 
set(handles.axes1,'OuterPosition',[0 28/h 1 (h-28)/h],'LooseInset',get(handles.axes1,'TightInset')+0.02)
% ========================================================================================


function pushbutton_Callback(hObject,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Save graph series in current axes
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname = uigetdir(cd,'Select directory to save image series.');
if pathname
    % Turn UI element visibility off
    figPos = get(handles.PlotFigure,'Position');
    h = figPos(4);
    set(handles.pushbutton1,'Visible','off');
    set(handles.slider1,'Visible','off');
    set(handles.axes1,'OuterPosition',[0 0 1 1])
    set(handles.PlotFigure,'Visible','off','PaperPositionMode','auto','InvertHardCopy', 'off')
    disp('Saving figure stack...')

    for i = 1:size(handles.Export.Measurement1,2)
        formplot(handles,i);
        print(handles.PlotFigure, [pathname,filesep,'scatter',numseq(i,3)], '-dpng',['-r',num2str(get(0, 'ScreenPixelsPerInch'))])
    end

    % Turn UI element visibility on
    set(handles.pushbutton1,'Visible','on');
    set(handles.slider1,'Visible','on');
    set(handles.axes1,'OuterPosition',[0 28/h 1 (h-28)/h]);
    set(handles.PlotFigure,'Visible','on')
    formplot(handles,get(handles.slider1,'Value'));
end
% ========================================================================================
