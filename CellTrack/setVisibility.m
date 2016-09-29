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
    set(handles.text2B,'String','N/A','ForegroundColor',[0.9 0.9 0.9]) % Used for display output
    set(handles.edit2B,'ForegroundColor',[0.9 0.9 0.9]) % Used for display output

else
    
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
        set(handles.text2B,'String','DIC','ForegroundColor',[13 84 171]/255)
        set(handles.edit2B,'ForegroundColor',[0 0 0]) % Used for display output

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
        set(handles.text2B,'String','Phase','ForegroundColor',[0 0 0])
    end
end


