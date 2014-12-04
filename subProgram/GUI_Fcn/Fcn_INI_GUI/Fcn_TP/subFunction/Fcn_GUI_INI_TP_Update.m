function Fcn_GUI_INI_TP_Update(varargin)
% this function is used to update ...
hObject = varargin{1};
handles = guidata(hObject);
global CI
%
% try
main = handles.MainGUI;                     % get the handle of OSCILOS_long
if(ishandle(main))
    mainHandles = guidata(main);            %
    % ------------------------------------
    String_Listbox=get(mainHandles.listbox_Info,'string');
    ind=find(ismember(String_Listbox,'<HTML><FONT color="blue">Information 2:'));
    nLength=size(String_Listbox);
    if isempty(ind)
        indStart=nLength(1);
    else
        indStart=ind-1;
        for i=nLength(1):-1:indStart+1 
            String_Listbox(i)=[];
        end
    end
    String_Listbox{indStart+1}='<HTML><FONT color="blue">Information 2:';
    String_Listbox{indStart+2}=['<HTML><FONT color="blue">Mean flow and thermal properties:'];
    % -----------------------------------
    switch CI.TP.isNoHA
        case 0
            set(mainHandles.INI_FM, 'Enable', 'on',...
                                    'visible','on');
            String_Listbox{indStart+3}=['The mean heat release rate is (are):'];
            String_Listbox{indStart+4}=[num2str(CI.TP.Q./1000) ' kW'];
        case 1
            set(mainHandles.INI_FM, 'Enable', 'off',...
                                    'visible','off');
            String_Listbox{indStart+3}=['There is no heat addition!'];
            String_Listbox{indStart+4}=[''];
    end
    String_Listbox{indStart+5}=['The mean temperature in different sections are:'];
    String_Listbox{indStart+6}=[num2str(CI.TP.T_mean(1,:)) ' K'];
    String_Listbox{indStart+7}=['The mean pressure in different sections are:'];
    String_Listbox{indStart+8}=[num2str(CI.TP.p_mean(1,:)) ' Pa'];
    String_Listbox{indStart+9}=['The mean velocity in different sections are:'];
    String_Listbox{indStart+10}=[num2str(CI.TP.u_mean(1,:)) ' m/s'];
    String_Listbox{indStart+11}=['The mean speed of sound in different sections are:'];
    String_Listbox{indStart+12}=[num2str(CI.TP.c_mean(1,:)) ' m/s'];
    set(mainHandles.listbox_Info,'string',String_Listbox,'value',1);
    guidata(hObject, handles);
else
end
% 
% ----------------------------end------------------------------------------