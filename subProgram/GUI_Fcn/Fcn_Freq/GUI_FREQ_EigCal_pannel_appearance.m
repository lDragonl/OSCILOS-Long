function GUI_FREQ_EigCal_pannel_appearance(varargin)
% This program is used to determine the appearance of GUI_FREQ_EigCal
%
% we need to check if this function has ever been run
% we still need to check which UI is visible or not, enable or disable
%
%
% 1. without heat perturbation
% 2. with heat perturbations
% 2.1 with linear heat perturbations
% 2.2 with nonlinear heat perturbations
% 2.2.1 with nonlinear flame model
% 2.2.2 flame describing functions are from experiment or CFD data
hObject = varargin{1};
handles = guidata(hObject);
global CI
switch CI.TP.isHP_final
    case 1    % with heat perturbations
        switch CI.indexFM
            case 0   % from flame model
                switch CI.FM.NL.style   % nonlinear index
                    case 1              % linear flame model

                        CI.EIG.APP_style = 12;
                    otherwise           % nonlinear flame model
                        CI.EIG.APP_style = 21;           % velocity ratios can be selected
                end
            case 1  % from exp or CFD flame describing function
                CI.EIG.APP_style = 22;                  % velocity ratios cannot be selected
        end
    case 0    % without perturbation
        CI.EIG.APP_style = 11;                   % without heat addition
end

%
switch CI.IsRun.GUI_FREQ_EigCal
    case 0
        set(handles.ObjEditEnable_AOC,    'Enable','off');
        switch CI.EIG.APP_style
            case {11,12}
                set(handles.ObjEditVisible_uRatio,      'visible','off');
            case 21                             % nonlinear flame model
                set(handles.ObjEditVisible_uRatio,      'visible','on');
            case 22                             % nonlinear flame model
                CI.EIG.FDF.uRatioSp = CI.FMEXP.uRatio;
                CI.EIG.FDF.uRatioNum = length(CI.EIG.FDF.uRatioSp);
                set(handles.edit_uRatio_min,       'string',num2str(CI.EIG.FDF.uRatioSp(1)));
                set(handles.edit_uRatio_max,       'string',num2str(CI.EIG.FDF.uRatioSp(end)));
                set(handles.edit_uRatio_SampNum,   'string',num2str(CI.EIG.FDF.uRatioNum));
        end
    case 1 
        set(handles.ObjEditEnable_AOC,    'Enable','on');
        GUI_FREQ_EigCal_PLOT(hObject)
        switch CI.EIG.APP_style
            case {11,12}
                set(handles.ObjEditVisible_uRatio,      'visible','off');
                set(handles.pop_PlotType,...
                'string',{'Map of eigenvalues';'Modeshape'},...
                'enable','on'); 
            case {21,22}                             % nonlinear flame model
                set(handles.ObjEditVisible_uRatio,      'visible','on');
                % -----------------------------------
                %  pannel related to uRatio 
                set(handles.slider_uRatio,...
                        'Enable','on',...
                        'min',1,...
                        'max',CI.EIG.FDF.uRatioNum,...
                        'value',1,...
                        'SliderStep',[1/(CI.EIG.FDF.uRatioNum-1), 1/(CI.EIG.FDF.uRatioNum-1)]);
                set(handles.edit_uRatio,'string',num2str(CI.EIG.FDF.uRatioSp(1)));
                set(handles.edit_uRatio_min,       'string',num2str(CI.EIG.FDF.uRatioSp(1)));
                set(handles.edit_uRatio_max,       'string',num2str(CI.EIG.FDF.uRatioSp(end)));
                set(handles.edit_uRatio_SampNum,   'string',num2str(CI.EIG.FDF.uRatioNum));
                set(handles.pop_PlotType,...
                'string',{'Map of eigenvalues';'Modeshape'; 'Evolution of eigenvalue with velocity ratio'},...
                'enable','on'); 
                % -----------------------------------
        end
        set(handles.pop_numMode,    'enable','on');  
        set(handles.pb_AOC_Plot,    'enable','on');
        eigenvalue = CI.EIG.Scan.EigValCol{1};
        for k = 1:length(eigenvalue)
           StringMode{k} = ['Mode number: ' num2str(k)]; 
        end
        set(handles.pop_numMode,    'string',StringMode)
        %
        data_num(:,1)=abs(imag(eigenvalue)./2./pi);
        data_num(:,2)=real(eigenvalue);
        data_cell=num2cell(data_num);
        set(handles.uitable,'data',data_cell);         % Update the table
        guidata(hObject, handles);
end
switch CI.EIG.APP_style
    case 22                             % from exp or CFD flame describing function
        set(handles.edit_uRatio_min,       'enable','off');
        set(handles.edit_uRatio_max,       'enable','off');
        set(handles.edit_uRatio_SampNum,   'enable','off');
end

%
% check the scan domain information
switch CI.IsRun.GUI_FREQ_EigCal_AD 
    case 0
        CI.EIG.Scan.FreqMin = 0;
        CI.EIG.Scan.FreqMax = 1000;
        CI.EIG.Scan.FreqNum = 10;
        CI.EIG.Scan.GRMin   = -200;
        CI.EIG.Scan.GRMax   = 200;
        CI.EIG.Scan.GRNum   = 10;
    case 1
end
guidata(hObject, handles);
%
assignin('base','CI',CI);


