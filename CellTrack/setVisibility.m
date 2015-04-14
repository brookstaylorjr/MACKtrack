function setVisibility(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Set GUI element visibility/text as appropriate, based on image type
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set visibility appropriately

if strcmp(handles.parameters.ImageType,'None')
    % Panel 6- Cell parameters. Drop entirely.
    set(handles.uipanel6,'Visible','off');
    set(handles.axes6A,'Visible','off')
    set(handles.pushbuttonSwitch2,'Visible','off') % Pushbutton "tab" for panel 6
    set(handles.text2B,'String','Aux') % Used for display output
    
    % Panel 5 - Nuclear parameters
    set(handles.slider5A,'Visible','off')
    set(handles.text5A,'Visible','off')
    set(handles.edit5A,'Visible','off')
    set(handles.text5B,'Visible','off')
    set(handles.edit5B,'Visible','off')
    axes_children = get(handles.axes5A,'children');
    set(axes_children,'Visible','off')
    set(handles.axes5A,'Visible','off')
else
    % Panel 5 - Nuclear parameters
    set(handles.slider5A,'Visible','on')
    set(handles.text5A,'Visible','on')
    set(handles.edit5A,'Visible','on')
    set(handles.text5B,'Visible','on')
    set(handles.edit5B,'Visible','on')
    set(handles.axes5A,'Visible','on')
    axes_children = get(handles.axes5A,'children');
    set(axes_children,'Visible','on')
    
    if strcmp(handles.parameters.ImageType,'DIC')
        % Panel 6 - Cell parameters. Set to DIC-specific fields
        set(handles.uipanel6_C,'Visible','off');
        set(handles.uipanel6_D,'Visible','on');
        set(handles.text6D,'Visible','off');
        set(handles.edit6D,'Visible','off');
        set(handles.text6D_2,'Visible','off');
        set(handles.uipanel6,'Title','DIC Parameters')
        set(handles.pushbuttonSwitch2,'Visible','on') % Pushbutton for panel 6
        set(handles.pushbuttonSwitch2,'String','DIC')
        set(handles.text2B,'String','DIC')
    else % (phase)
        % Panel 6- Cell parameters. Set to phase-specific fields
        set(handles.uipanel6_C,'Visible','on');
        set(handles.uipanel6_D,'Visible','off');
        set(handles.text6D,'Visible','on');
        set(handles.edit6D,'Visible','on');
        set(handles.text6D_2,'Visible','on');
        set(handles.uipanel6,'Title','Phase Parameters')
        set(handles.pushbuttonSwitch2,'Visible','on') % Pushbutton for panel 6
        set(handles.pushbuttonSwitch2,'String','Phase')
        set(handles.text2B,'String','Phase')
    end
end


