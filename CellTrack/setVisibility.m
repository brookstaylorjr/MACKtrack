function setVisibility(handles)
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
% Set GUI element visibility/text as appropriate, based on image type
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set visibility appropriately


switch lower(handles.parameters.ImageType)

    case 'none'  % Drop panel 6 (cell parameters) entirely.
        set(handles.uipanel6,'Visible','off');
        set(handles.axes6A,'Visible','off')
        set(handles.pushbuttonSwitch2,'Visible','off') % Pushbutton "tab" for panel 6
        set(handles.text2B,'String','N/A','ForegroundColor',[0.9 0.9 0.9]) % Used for display output
        set(handles.edit2B,'ForegroundColor',[0.9 0.9 0.9]) % Used for display output
        set(handles.text2I,'ForegroundColor',handles.gray);
    
    case 'dic' % Set to DIC-specific fields
        set(handles.axes6A,'Visible','on');
        set(handles.uipanel6_A,'Visible','on');
        set(handles.uipanel6_E,'Visible','off');

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
        set(handles.text2I,'ForegroundColor','r');


    case 'phase' % Set to phase-contrast  fields
        set(handles.axes6A,'Visible','on');
        set(handles.uipanel6_A,'Visible','on');
        set(handles.uipanel6_E,'Visible','off');
        
        set(handles.uipanel6_C,'Visible','on');
        set(handles.uipanel6_D,'Visible','off');
        set(handles.text6D,'Visible','on');
        set(handles.edit6D,'Visible','on');
        set(handles.text6D_2,'Visible','on');
        set(handles.uipanel6,'Title','Phase Parameters')
        set(handles.pushbuttonSwitch2,'Visible','on') % Pushbutton for panel 6
        set(handles.pushbuttonSwitch2,'String','Phase')
        set(handles.text2B,'String','Phase','ForegroundColor',[0 0 0])
        set(handles.text2I,'ForegroundColor','r');
        
    case 'fluorescence' % Set to fluorescence-specific fields
        set(handles.axes6A,'Visible','off');
        set(handles.uipanel6_A,'Visible','off');
        set(handles.uipanel6_E,'Visible','on');
        
        set(handles.uipanel6_C,'Visible','off');
        set(handles.uipanel6_D,'Visible','on');
        set(handles.text6D,'Visible','off');
        set(handles.edit6D,'Visible','off');
        set(handles.text6D_2,'Visible','off');
        set(handles.uipanel6,'Title','Cell Parameters')
        set(handles.pushbuttonSwitch2,'Visible','on') % Pushbutton for panel 6
        set(handles.pushbuttonSwitch2,'String','Cell')
        set(handles.text2B,'String','Cell','ForegroundColor',[13 84 171]/255)
        set(handles.edit2B,'ForegroundColor',[0 0 0]) % Used for display output
        set(handles.text2I,'ForegroundColor','r');
        

end


