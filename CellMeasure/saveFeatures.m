function featuresOut = saveFeatures(imageIn,featuresIn)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% SAVEFEATURES allows user selection of features in a given image for use in featureModule
% Features are saved from images in the featuresOut structure, named feature001, 
% feature002, etc.
%
% imageIn         input image used in feature selection
% featuresIn      (optional) input cell array of all saved subimages
%
% featuresOut     output cell array of all saved subimages
%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


%  Input checking: create cell matrix if one is not provided.
if nargin<2
    featuresIn = cell(0);
end

% Reads (opens) input image and assigns it to a variable, "origImage"
origImage = double(imread(imageIn));
origImage = (origImage- min(origImage(:)))/(max(origImage(:))- min(origImage(:))); % Scale image 0 to 1

% Display the image in a new figure window.
h = figure;
imshow(origImage,'InitialMagnification','fit'), colormap jet
findall(h,'Type','uitoolbar')

title('Select Upper left and lower-right corner of desired feature(s)')
hAx = get(h,'CurrentAxes');
orig_axis = axis(hAx);
hZoom = zoom(gcf);
% Allow user to specify supregions of image.
loopControl = 1;

while loopControl == 1
    %%% Manual outline of the object, see local function 'getpoints' below.
    %%% Continue until user has drawn a valid shape
    [x,y] = getpoints(hAx,hZoom,orig_axis);
    
    %nrows,ncols] = size(origImage);
    %[X,Y] = meshgrid(1:ncols,1:nrows);
    %objMask = inpolygon(X,Y, x,y);
    %outputLabel(objMask) = i;
    featuresIn = cat(1,featuresIn,origImage(min(y):max(y),min(x):max(x)));
    ButtonName=questdlg('Would you like to outline another feature in this image?', ...
        'Manual Trace', ...
        'Yes', 'No', 'Yes');
    if(strcmp(ButtonName, 'No'))
        loopControl = 0;
    end
end
close(h)
featuresOut = featuresIn;


% ========================================================================================

function [x,y] = getpoints(AxisHandle,ZoomHandle, OriginalAxis)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% GETPOINTS allow user to select upper left/lower right of feature- modeled closely after
% 'getpoints' subfunction in CellProfiler
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

Position = get(AxisHandle,'Position');
FigureHandle = (get(AxisHandle, 'Parent'));
PointHandles = [];
x = [];
y = [];
NbrOfPoints = 0;
done = 0;

ImageHandle = get(AxisHandle,'children');
set(ImageHandle,'ButtonDownFcn','');
hold on
while ~done;
    UserInput = waitforbuttonpress;                            % Wait for user input
    SelectionType = get(FigureHandle,'SelectionType');         % Get information about the last button press
    CharacterType = get(FigureHandle,'CurrentCharacter');      % Get information about the character entered
    ZoomState = get(ZoomHandle,'Enable');
    
    % Process zoom state vs selection state
    if strcmp(ZoomState,'off') 
        % 1) Normal left mouse button + 0 or 1 points: add new point.
        if UserInput == 0 && strcmp(SelectionType,'normal') && NbrOfPoints < 2
            % Get the new point and store it
            CurrentPoint  = get(AxisHandle, 'CurrentPoint');
            x = [x round(CurrentPoint(2,1))];
            y = [y round(CurrentPoint(2,2))];
            NbrOfPoints = NbrOfPoints + 1;

            % Plot the new point
            h = plot(round(CurrentPoint(2,1)),round(CurrentPoint(2,2)),'r.');
            set(AxisHandle,'Position',Position)                   %  Keep Matlab from moving title text/resizing image when first point is plotted
            PointHandles = [PointHandles h];
        
        % 2) Right mousebutton/backspace key + 1 or 2 points: delete last point
        elseif NbrOfPoints > 0 && ((UserInput == 0 && strcmp(SelectionType,'alt')) || (UserInput == 1 && CharacterType == char(8))) % (ASCII code for backspace is 8)
            NbrOfPoints = NbrOfPoints - 1;
            x = x(1:end-1);
            y = y(1:end-1);
            delete(PointHandles(end));
            PointHandles = PointHandles(1:end-1);

        % 3) Enter: save object and draw 
        elseif (NbrOfPoints == 2) && ((UserInput == 1 && CharacterType == char(13)) || (UserInput == 0 && strcmp(SelectionType,'open'))) && strcmp(ZoomState,'off')
            % Indicate that we are done
            done = 1;
            % Remove plotted points
            if ~isempty(PointHandles)
                delete(PointHandles)
            end
        end
    
        
        else % - - - - - - Zoom on - - - - - - - - -
         
         % 1) Left mouse button: zoom in/out by factor of 2
        if UserInput == 0 && strcmp(SelectionType,'normal')    
            CurrentPoint  = get(AxisHandle, 'CurrentPoint');
            xzoom = round(CurrentPoint(2,1));
            yzoom = round(CurrentPoint(2,2));
            if strcmp(get(ZoomHandle,'Direction'),'in')
                zoomcenter(xzoom,yzoom,2);
            else
                zoomcenter(xzoom,yzoom,0.5);
            end
        
        % 5) Double click or enter: reset original axes
        elseif ((UserInput == 1 && CharacterType == char(13)) || (UserInput == 0 && strcmp(SelectionType,'open')))
            axis(gca,OriginalAxis);
    
        end
        
    end

    
    % Remove old box variable and draw new
    if exist('featurebox','var')
        delete(featurebox) % Delete the graphics object
        clear featurebox % Clear the variable
    end
    
    if NbrOfPoints > 1
        total_xpts = max(x)-min(x)+1;
        total_ypts = max(y)-min(y)+1;
        xpts = [min(x)*ones(1,total_ypts), min(x):max(x), max(x)*ones(1,total_ypts),min(x):max(x)];
        ypts = [min(y):max(y), min(y)*ones(1,total_xpts),min(y):max(y),max(y)*ones(1,total_xpts)];
        featurebox = plot(xpts,ypts,'r');
        drawnow
    end
end
hold off
% ========================================================================================


function zoomcenter(x,y,factor)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%ZOOMCENTER Zoom in and out of a specified point on a 2-D plot.
%
% x       x position of zoom center
% y       y position of zoom center
% factor  zoom factor (identical to MATLAB's zoom)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ax = gca;
cax = axis(ax);
daxX = (cax(2)-cax(1))/factor(1)/2;
daxY = (cax(4)-cax(3))/factor(end)/2;
axis(ax,[x+[-1 1]*daxX y+[-1 1]*daxY]);
