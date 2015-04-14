function histogramSeries(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Make histogram for frame- initialize "cycle frame" slider
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Initialize axes and slider and set properties
if isempty(get(handles.PlotFigure,'Children'))
    % Axes1: measurement 1
    handles.axes1 = axes('Parent',handles.PlotFigure);
    % Axes1: measurement 2
    handles.axes2 = axes('Parent',handles.PlotFigure);
    % Initialize slider
    stepSize = 1/size(handles.Export.Measurement1,2);
    handles.slider1  = uicontrol('Style', 'slider','Max',size(handles.Export.Measurement1,2),...
        'Min',1,'Value',1,'SliderStep', [stepSize stepSize*10]);
    set(handles.slider1,'Callback',{@slider_Callback,handles},'BackgroundColor',[.9 .9 .9])
    % Initialize "save" button (single)
    handles.pushbutton1  = uicontrol('Style', 'pushbutton','String','Save Top');
    % Initialize "save" button (double)
    handles.pushbutton2  = uicontrol('Style', 'pushbutton','String','Save Both');
    % Callback/resize functions
    set(handles.pushbutton1,'Callback',{@pushbutton1_Callback,handles},'BackgroundColor',[.9 .9 .9])
    set(handles.pushbutton2,'Callback',{@pushbutton2_Callback,handles},'BackgroundColor',[.9 .9 .9])
    set(handles.PlotFigure,'ResizeFcn',{@fig_resize,handles})
    
end

frame = get(handles.slider1,'Value');

formplot(handles,frame,1)
formplot(handles,frame,2)

% ========================================================================================

function formplot(handles, frame, axesnum)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Form histogram in selected axes
%
% handles     master structure
% frame       current time point
% axesnum     target axes for graphs (1 or 2)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% Update time element
pos = handles.Export.Measurement1Info.Times(frame);
timeHr = fix(pos);
timeMin = round(60*(pos-timeHr));

% Construct field names/quantities/titles for current axes
measurementinfo = ['Measurement',num2str(axesnum),'Info'];
measurement = ['Measurement',num2str(axesnum)];
axesname = ['axes',num2str(axesnum)];
num_bins = handles.Options.(['HistogramBins',num2str(axesnum)]);
if axesnum==1
    titlestr = [numseq(abs(timeHr),2),':',numseq(abs(timeMin),2)];
        if (timeMin<0) || (timeHr<0)
            titlestr = ['-',titlestr];
        end
    barcolor = [194,224,120]/255; % light blue-gray
else
    titlestr = ' ';
    barcolor = [84,123,151]/255; % gray-blue
end

data_range = handles.Export.(measurementinfo).MeasurementBounds(1):...
    diff(handles.Export.(measurementinfo).MeasurementBounds)/num_bins:handles.Export.(measurementinfo).MeasurementBounds(2);

% Make histogram; convert to relative frequency
n = hist(handles.Export.(measurement)(:,frame),data_range);
n = n./sum(n);
bar(handles.(axesname),data_range,n,'BarWidth',0.88,'EdgeColor','none','FaceColor',barcolor)
axis(handles.(axesname),[min(data_range),max(data_range),0, handles.Options.(['HistogramMax',num2str(axesnum)])])


% Set title, axis labels and properties
title(handles.(axesname),titlestr,'FontSize',16)
xlabel(handles.(axesname),[handles.Export.(measurementinfo).Name],'FontSize',14);
ylabel(handles.(axesname),'Relative Frequency','FontSize',14);
set(handles.(axesname),'XTick',handles.Export.(measurementinfo).MeasurementTicks,'XTickLabel',handles.Export.(measurementinfo).MeasurementTickLabels,...
'YTick',0.1:0.1:0.5,'Box','off','TickLength',[0.005 0.005]);


% Make overlay
% Gaussian (1-D)
if (strcmp(handles.Options.Grouping,get(handles.radiobutton8C,'String')) && (axesnum==(handles.Options.Measure2Grouping+1)))
    hold(handles.(axesname),'on')
    % Solid line: PDF of combined function
    pdf_x = (min(data_range) : (max(data_range)-min(data_range))/999: max(data_range))';
    pdf_y = pdf(handles.Export.Gaussian.Dist{frame},pdf_x);
    stepmult = (data_range(2)-data_range(1))/(pdf_x(2)-pdf_x(1));
    pdf_y = pdf_y/sum(pdf_y)*stepmult;
    plot(handles.(axesname),pdf_x,pdf_y,'Color',[0.2 0.2 0.2],'LineWidth',2);
    % Dashed lines: PDFs of each subpopulation
    pdf_subpop = cell(1,handles.Options.Subpopulations);
    for i = 1:handles.Options.Subpopulations
        pdf_subpop{i} = normpdf(pdf_x,handles.Export.Gaussian.Dist{frame}.mu(i),sqrt(handles.Export.Gaussian.Dist{frame}.Sigma(:,:,i)));
        pdf_subpop{i} = pdf_subpop{i}/sum(pdf_subpop{i})*stepmult*handles.Export.Gaussian.Dist{frame}.PComponents(i);
        plot(handles.(axesname),pdf_x,pdf_subpop{i},'--','Color',[0.2 0.2 0.2],'LineWidth',2);
    end    
    hold(handles.(axesname),'off')
% 2D Gaussian case
elseif strcmp(handles.Options.Grouping,get(handles.radiobutton8E,'String'))    
    hold(handles.(axesname),'on')
    % Solid line: PDF of combined function
    pdf_x = (min(data_range) : (max(data_range)-min(data_range))/999: max(data_range))';
    %pdf_y = pdf(handles.Export.Gaussian.Dist{frame},pdf_x);
    stepmult = (data_range(2)-data_range(1))/(pdf_x(2)-pdf_x(1));
    %pdf_y = pdf_y/sum(pdf_y)*stepmult;
    
    % Dashed lines: PDFs of each subpopulation
    pdf_subpop = cell(1,handles.Options.Subpopulations);
    for i = 1:handles.Options.Subpopulations
        pdf_subpop{i} = normpdf(pdf_x,handles.Export.Gaussian.Dist{frame}.mu(i,axesnum),sqrt(handles.Export.Gaussian.Dist{frame}.Sigma(axesnum,axesnum,i)));
        pdf_subpop{i} = pdf_subpop{i}/sum(pdf_subpop{i})*stepmult*handles.Export.Gaussian.Dist{frame}.PComponents(i);
        plot(handles.(axesname),pdf_x,pdf_subpop{i},'--','Color',[0.2 0.2 0.2],'LineWidth',2);
    end    
    hold(handles.(axesname),'off')

end




% ========================================================================================



function slider_Callback(hObject,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Update time with slider movement/ redraw graph
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
set(hObject,'Value',round(get(hObject,'Value')))
histogramSeries(handles)
% ========================================================================================

function fig_resize(~,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Resize histograms and slider smoothly
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
figPos = get(handles.PlotFigure,'Position');
h = figPos(4);
h_adj = h-28;
w = figPos(3);
set(handles.axes1,'OuterPosition',[0 (28+h_adj/2)/h 1 (h_adj/2)/h],'LooseInset',get(handles.axes1,'TightInset')+0.02)
set(handles.axes2,'OuterPosition',[0 28/h 1 (h_adj/2)/h],'LooseInset',get(handles.axes2,'TightInset')+0.02)
set(handles.slider1,'Position',[w/4-70 5 w/2 20]);
set(handles.pushbutton1,'Position',[3*w/4-60 4 80 24]); 
set(handles.pushbutton2,'Position',[3*w/4+25 4 80 24]); 
% ========================================================================================

function pushbutton1_Callback(hObject,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Save graph series in current axes (top histogram only)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname = uigetdir(cd,'Select directory to save image series.');

if pathname
    disp('Saving figure stack...')
    % Create temporary figure and move axes1 over
    figPos = get(handles.PlotFigure,'Position');
    h = figPos(4);
    h_adj = h-28;
    handles.tmpfig = figure;
    set(handles.tmpfig,'DockControls','off','Color',[1 1 1],'Visible','off','PaperPositionMode','auto','InvertHardCopy', 'off','Position',figPos./[1 1 1 2])
    set(handles.axes1,'Parent',handles.tmpfig,'OuterPosition',[0 0 1 1],'LooseInset',get(handles.axes1,'TightInset')+0.02)

    % Cycle through, updating axes and saving
    for i = 1%:size(handles.Export.Measurement1,2)
        formplot(handles,i,1)
        print(handles.tmpfig, [pathname,filesep,'hist',numseq(i,3)], '-dpng',['-r',num2str(get(0, 'ScreenPixelsPerInch'))])
    end

    % Move axes back to PlotFigure, reset size, and close tmpfig
    set(handles.axes1,'Parent',handles.PlotFigure);
    set(handles.axes1,'OuterPosition',[0 (28+h_adj/2)/h 1 (h_adj/2)/h],'LooseInset',get(handles.axes1,'TightInset')+0.02)
    close(handles.tmpfig)
    formplot(handles,get(handles.slider1,'Value'),1);
    disp('Done.')
end
% ========================================================================================



function pushbutton2_Callback(hObject,~,handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Save graph series in current axes (both histograms)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname = uigetdir(cd,'Select directory to save image series.');
if pathname
    disp('Saving figure stack...')
    % Turn UI element visibility off
    figPos = get(handles.PlotFigure,'Position');
    h = figPos(4);
    h_adj = h-28;
    set(handles.pushbutton1,'Visible','off');
    set(handles.pushbutton2,'Visible','off');
    set(handles.slider1,'Visible','off');
    set(handles.axes1,'OuterPosition',[0 0.5 1 0.5],'LooseInset',get(handles.axes1,'TightInset')+0.02)
    set(handles.axes2,'OuterPosition',[0 0 1 0.5],'LooseInset',get(handles.axes2,'TightInset')+0.02)
    set(handles.PlotFigure,'Visible','off','PaperPositionMode','auto','InvertHardCopy', 'off')
    

    for i = 1:size(handles.Export.Measurement1,2)
        formplot(handles,i,1)
        formplot(handles,i,2)
        print(handles.PlotFigure, [pathname,filesep,'double_hist',numseq(i,3)], '-dpng',['-r',num2str(get(0, 'ScreenPixelsPerInch'))])
    end

    % Turn UI element visibility back on
    set(handles.pushbutton1,'Visible','on');
    set(handles.pushbutton2,'Visible','on');
    set(handles.slider1,'Visible','on');
    set(handles.axes1,'OuterPosition',[0 (28+h_adj/2)/h 1 (h_adj/2)/h],'LooseInset',get(handles.axes1,'TightInset')+0.02)
    set(handles.axes2,'OuterPosition',[0 28/h 1 (h_adj/2)/h],'LooseInset',get(handles.axes2,'TightInset')+0.02)
    set(handles.PlotFigure,'Visible','on')
    formplot(handles,get(handles.slider1,'Value'),1);
    formplot(handles,get(handles.slider1,'Value'),2);
    disp('Done.')
end
% ========================================================================================
