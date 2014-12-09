function Fcn_GUI_INI_TP_Update(varargin)
% this function is used to update ...
hObject = varargin{1};
handles = guidata(hObject);
global CI
% -------------identify if the GUI_INI_FM.m need to be lauched ------------
switch CI.CD.isHA
    case 1   % With heat addition
        CI.TP.isHA_final = 1;               % With heat addition
        CI.TP.isHP_final = 1;               % With heat perturbations
        % ------Exceptions---------
        % case 1: without heat perturbations
        if CI.CD.isHP == 0   %  in case without heat perturbations
            CI.TP.isHP_final = 0; 
        end
        % case 2: all TRatio == 1, and consider there is no heat addition
        % and perturbations
        if CI.TP.TRatio == 1
            CI.TP.isHA_final = 0;
            CI.TP.isHP_final = 0;
        end
    case 0   % Without heat addition
        CI.TP.isHA_final = 0;               % No heat addition
        CI.TP.isHP_final = 0;               % No heat perturbation
end
assignin('base','CI',CI);
%
% ---------------------Update other GUI window-----------------------------
%
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
    if CI.TP.isHA_final == 0    % No heat addition
        String_Listbox{indStart+3}=['There is no heat addition!'];
        String_Listbox{indStart+4}=[''];
    elseif  CI.TP.isHA_final == 0 
        String_Listbox{indStart+3}=['The mean heat release rate is (are):'];
        String_Listbox{indStart+4}=[num2str(CI.TP.Q./1000) ' kW'];
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
    % -----------
    if CI.TP.isHP_final == 0    % No heat perturbations
        set(mainHandles.INI_FM, 'Enable', 'off',...
                                'visible','off');
        set(mainHandles.INI_BC, 'Enable', 'on');
        Fcn_default_flame_model   
    elseif CI.TP.isHP_final == 1 % with heat perturbations
        set(mainHandles.INI_FM, 'Enable', 'on',...
                                'visible','on');
        set(mainHandles.INI_BC, 'Enable', 'off');
    end
    guidata(hObject, handles);
else
end

function Fcn_default_flame_model
global CI
% when the pannel for setting flame model is not avaliable, this function
% is launched to initialize the parameters corresponding to the flame model
% 
CI.FM.FTF.num = 0;
CI.FM.FTF.den = 1;
CI.FM.FTF.tauf = 0;
assignin('base','CI',CI);
% ----------------------------end------------------------------------------