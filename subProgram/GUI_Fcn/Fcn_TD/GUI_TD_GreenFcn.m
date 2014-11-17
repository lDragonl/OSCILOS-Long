function varargout = GUI_TD_GreenFcn(varargin)
% GUI_TD_GreenFcn MATLAB code for GUI_TD_GreenFcn.fig
%      GUI_TD_GreenFcn, by itself, creates a new GUI_TD_GreenFcn or raises the existing
%      singleton*.
%
%      H = GUI_TD_GreenFcn returns the handle to a new GUI_TD_GreenFcn or the handle to
%      the existing singleton*.
%
%      GUI_TD_GreenFcn('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI_TD_GreenFcn.M with the given input arguments.
%
%      GUI_TD_GreenFcn('Property','Value',...) creates a new GUI_TD_GreenFcn or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before GUI_TD_GreenFcn_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to GUI_TD_GreenFcn_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help GUI_TD_GreenFcn

% Last Modified by GUIDE v2.5 24-Oct-2014 09:17:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @GUI_TD_GreenFcn_OpeningFcn, ...
                   'gui_OutputFcn',  @GUI_TD_GreenFcn_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
%
% -------------------------------------------------------------------------
%          
% --- Executes just before GUI_TD_GreenFcn is made visible.
function GUI_TD_GreenFcn_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to GUI_TD_GreenFcn (see VARARGIN)
indexEdit = 0;
switch indexEdit 
    case 0
        %--------------------------------------------------------------------------
        dontOpen = false;
        mainGuiInput = find(strcmp(varargin, 'OSCILOS_long'));
        if (isempty(mainGuiInput)) ...
            || (length(varargin) <= mainGuiInput) ...
            || (~ishandle(varargin{mainGuiInput+1}))
            dontOpen = true;
        else % load from the main GUI
            % handles of main GUI
            handles.MainGUI = varargin{mainGuiInput+1};
            try
                handles.ExampleGUI = varargin{mainGuiInput+2};
            catch
            end
            % Obtain handles using GUIDATA with the caller's handle 
            mainHandles         = guidata(handles.MainGUI);
            % background colors
            handles.bgcolor     = mainHandles.bgcolor;
            % fontsize
            handles.FontSize    = mainHandles.FontSize;
            %
            handles.sW          = mainHandles.sW;
            handles.sH          = mainHandles.sH;
            handles.indexApp    = 0;
            % Update handles structure
            guidata(hObject, handles);
            % Initialization
            GUI_Pannel_Initialization(hObject, eventdata, handles)
        end
        % Update handles structure
        guidata(hObject, handles);
        if dontOpen
           disp('-----------------------------------------------------');
           disp('This is a subprogram. It cannot be run independently.') 
           disp('Please load the program "OSCILOS_long'' from the ')
           disp('parent directory!')
           disp('-----------------------------------------------------');
        else
           uiwait(hObject);
        end
    case 1
        global CI
        %  only used for debugging the code
        handles.bgcolor{1} = [0.95, 0.95, 0.95];
        handles.bgcolor{2} = [0, 0, 0];
        handles.bgcolor{3} = [.75, .75, .75];
        handles.bgcolor{4} = [0.90,0.90,1];
        %
        handles.sW  = 800;
        handles.sH  = 600;
        %
        if ispc
            handles.FontSize(1)=11;                 % set the default fontsize
            handles.FontSize(2)=9;
        else
            handles.FontSize(1)=12;                 % set the default fontsize
            handles.FontSize(2)=10;   
        end
        CI.IsRun.GUI_TD_Convg = 0;
        CI.BC.num1      = 100;
        CI.BC.den1      = [1,100];
        CI.BC.num2      = 1;
        CI.BC.den2      = 1;
        CI.FM.FTF.num   = 200;
        CI.FM.FTF.den   = [1,200];
        assignin('base','CI',CI);                   % save the current information to the works
        handles.indexApp = 0;
        guidata(hObject, handles);  
        GUI_Pannel_Initialization(hObject, eventdata, handles)
        uiwait(hObject);
end
%
% -------------------------------------------------------------------------
%
function GUI_Pannel_Initialization(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI        
assignin('base','CI',CI);                   % save the current information to the workspace
%
set(0, 'units', 'points');
screenSize  = get(0, 'ScreenSize');                             % get the screen size
sW          = handles.sW;                                       % screen width
sH          = handles.sH;                                       % screen height
FigW        = sW.*1/2;                                          % window width
FigH        = sH.*3/5;                                          % window height
set(handles.figure,     'units', 'points',...
                        'position',[(screenSize(3)-FigW)./2 (screenSize(4)-FigH)./2 FigW FigH],...
                        'name','Green''function examination',...
                        'color',handles.bgcolor{3});
%----------------------------------------
% pannel axes
set(handles.uipanel_axes,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*8/20 FigW*19/20 FigH*11.5/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3});  
pannelsize=get(handles.uipanel_axes,'position');
pW=pannelsize(3);
pH=pannelsize(4);                
set(handles.axes1,      'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*1.5/10 pH*2.0/10 pW*7.5/10 pH*6/10],...
                        'fontsize',handles.FontSize(1),...
                        'color',handles.bgcolor{1},...
                        'box','on');                     
guidata(hObject, handles);
set(handles.pop_type,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*9/10 pW*5.5/10 pH*0.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',{  'Inlet boundary condition';...
                                    'Outlet boundary condition';...
                                    'Flame transfer function'},...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','left',...
                        'enable','on',...
                        'value',1);  
%----------------------------------------
% pannel input
set(handles.uipanel_Input,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*2.5/20 FigW*19/20 FigH*5/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_Input,'position');
pW=pannelsize(3);
pH=pannelsize(4);                       
%
set(handles.text_a1,    'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*6/10 pW*5/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Suggested saturation time [ms]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_a1,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*6/10 pW*3/10 pH*2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',1,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on',...
                        'enable','off');  
set(handles.text_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.5/10 pH*2/10 pW*5/10 pH*1.5/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','User defined saturation time [ms]:',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'horizontalalignment','left',...
                        'visible','on');                   
set(handles.edit_a2,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*6/10 pH*2/10 pW*3/10 pH*2/10],...
                        'fontsize',handles.FontSize(2),...
                        'string',0,...
                        'backgroundcolor',handles.bgcolor{1},...
                        'horizontalalignment','right',...
                        'visible','on');  

%----------------------------------------
% pannel AOC                   
set(handles.uipanel_AOC,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[FigW*0.5/20 FigH*0/20 FigW*19/20 FigH*2/20],...
                        'Title','',...
                        'visible','on',...
                        'highlightcolor',handles.bgcolor{3},...
                        'borderwidth',1,...
                        'fontsize',handles.FontSize(2),...
                        'backgroundcolor',handles.bgcolor{3}); 
pannelsize=get(handles.uipanel_AOC,'position');                    
pW=pannelsize(3);
pH=pannelsize(4);                
set(handles.pb_SaveFig,...
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*0.75/10 pH*2/10 pW*2.5/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Save figure',...
                        'backgroundcolor',handles.bgcolor{3},...
                        'enable','on'); 
set(handles.pb_OK,      'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*4/10 pH*2/10 pW*2.5/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','OK',...
                        'backgroundcolor',handles.bgcolor{3});
set(handles.pb_Cancel,....
                        'units', 'points',...
                        'Fontunits','points',...
                        'position',[pW*7.25/10 pH*2/10 pW*2.5/10 pH*6/10],...
                        'fontsize',handles.FontSize(2),...
                        'string','Cancel',...
                        'backgroundcolor',handles.bgcolor{3});
%
guidata(hObject, handles);  
Fcn_Pre_calculation(hObject) 
handles = guidata(hObject);
guidata(hObject, handles);  
Fcn_Update_Plots(hObject);
handles = guidata(hObject);
guidata(hObject, handles); 
%
% -------------------------------------------------------------------------
%
function Fcn_Pre_calculation(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global CI
global Green
switch CI.IsRun.GUI_TD_Convg
    case 0
        Green.num{1} = CI.BC.num1;
        Green.den{1} = CI.BC.den1;
        Green.num{2} = CI.BC.num2;
        Green.den{2} = CI.BC.den2;
        Green.num{3} = CI.FM.FTF.num;
        Green.den{3} = CI.FM.FTF.den;
        for ss = 1:3
            num     = Green.num{ss};
            den     = Green.den{ss};
            Nnum    = length(num);
            Nden    = length(den);
            if Nnum <=1 && Nden<=1
                Green.indexConst(ss)   = 1;
                Green.y1{ss}           = num./den;
                Green.t1{ss}           = 0;
                Green.tauConv1(ss)     = 0;
                Green.y2{ss}           = num./den;
                Green.t2{ss}           = 0;
                Green.tauConv2(ss)     = 0;
            else
                Green.indexConst(ss)   = 0;
                Green.sys{ss}          = tf(num,den);
                [Green.y1{ss},Green.t1{ss}]   = impulse(Green.sys{ss});
                t1                          = Green.t1{ss};
                Green.tauConv1(ss)     = t1(end);
                [Green.y2{ss},Green.t2{ss}]   = impulse(Green.sys{ss});
                t2                          = Green.t2{ss};
                Green.tauConv2(ss)     = t2(end);
            end
        end
        CI.TD.Green    = Green;
    case 1
        Green  = CI.TD.Green;
end
set(handles.pop_type, 'value',1)
set(handles.edit_a1, 'string', num2str(Green.tauConv1(1).*1e3));
set(handles.edit_a2, 'string', num2str(Green.tauConv2(1).*1e3));
assignin('base','CI',CI);
assignin('base','Green',Green);
guidata(hObject, handles) 
%
% -------------------------------------------------------------------------
%
function Fcn_Update_Plots(varargin)
hObject = varargin{1};
handles = guidata(hObject);
global Green
ss = get(handles.pop_type,'Value');
switch Green.indexConst(ss)
    case 0
        set(handles.edit_a2,    'enable','on');
        Green.tauConv2(ss) = str2double(get(handles.edit_a2, 'string'))./1e3;
        [Green.y2{ss},Green.t2{ss}]   = impulse(Green.sys{ss},Green.tauConv2(ss));
    case 1
        set(handles.edit_a2,    'enable','off');
            string={...
            'The transfer function is a constant value!';...
            'Examination of the Green''s function is not necessary!'};
        helpdlg(string,'')
end
assignin('base','Green',Green);
guidata(hObject, handles)
Fcn_Plots(hObject)
%
% -------------------------------------------------------------------------
%
function Fcn_Plots(varargin)
hObject     = varargin{1};
handles     = guidata(hObject);
hAxes1      = handles.axes1;
fontSize1   = handles.FontSize(1);
ss          = get(handles.pop_type,'Value');
global Green
try
    delete(hl)
catch
end
cla(hAxes1)
axes(hAxes1)
hold on
plot(hAxes1,Green.t1{ss}.*1e3,Green.y1{ss},'-','color','b','Linewidth',2); 
plot(hAxes1,Green.t2{ss}.*1e3,Green.y2{ss},'--','color','r','Linewidth',1); 
set(hAxes1,'YColor','k','Box','on','ygrid','on','xgrid','on');
set(hAxes1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
xlabel(hAxes1,'Time [ms]','Color','k','Interpreter','LaTex','FontSize',fontSize1);
ylabel(hAxes1,'Time domain response [-]','Color','k','Interpreter','LaTex','FontSize',fontSize1)
hold off
hl = legend(hAxes1,'Suggested','User defined');
set(hl,'interpreter','latex','Fontsize',handles.FontSize(2),'box','off','Unit','points')
%
% -------------------------------------------------------------------------
%
function pop_type_Callback(hObject, eventdata, handles)
% hObject    handle to pop_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Green
ss = get(handles.pop_type,'Value');
set(handles.edit_a1, 'string', num2str(Green.tauConv1(ss).*1e3));
set(handles.edit_a2, 'string', num2str(Green.tauConv2(ss).*1e3));
Fcn_Update_Plots(hObject)
handles = guidata(hObject);
guidata(hObject, handles);
%
% -------------------------------------------------------------------------
%
% --- Outputs from this function are returned to the command line.
function varargout = GUI_TD_GreenFcn_OutputFcn(hObject, eventdata, handles) 
varargout{1} = [];
delete(hObject);
%
% -------------------------------------------------------------------------
%
% --- Executes on button press in pb_OK.
function pb_OK_Callback(hObject, eventdata, handles)
% hObject    handle to pb_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CI
global Green
CI.TD.Green                = Green;
CI.IsRun.GUI_TD_Convg   = 1;
assignin('base','CI',CI);
uiresume(handles.figure);
%
% -------------------------------------------------------------------------
%
% --- Executes on button press in pb_Cancel.
function pb_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pb_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure);
%
% -------------------------------------------------------------------------
%
% --- Executes on button press in pb_SaveFig.
function pb_SaveFig_Callback(hObject, eventdata, handles)
handles     = guidata(hObject);
Fig         = figure;
copyobj(handles.axes1, Fig);
set(Fig,        'units','points')
posFig      = get(handles.figure,'position');
hAxes       = get(Fig,'children');
set(hAxes(1),       'units','points',...
                    'position',[60 60 200 150],...
                    'ActivePositionProperty','position')
posAxesOuter = [0 0 300 300];
set(Fig,        'units','points',...
                'position', [posFig(1)+0.5*posFig(3)-0.5*posAxesOuter(3),...
                            posFig(2)+0.5*posFig(4)-0.5*posAxesOuter(4),...
                            posAxesOuter(3:4)]) 
hl=legend(hAxes(1),'Suggested','User defined');
set(hl,'interpreter','latex','Fontsize',handles.FontSize(2),'box','off','Unit','points')
%
% -------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function pop_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
% -------------------------------------------------------------------------
%
function edit_a1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datEdit = str2double(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
%
% -------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_a1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
% -------------------------------------------------------------------------
%
function edit_a2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datEdit = str2double(get(hObject, 'String'));
if isnan(datEdit)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
Fcn_Update_Plots(hObject, eventdata, handles)
%
% -------------------------------------------------------------------------
%
% --- Executes during object creation, after setting all properties.
function edit_a2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
% -------------------------------------------------------------------------
%
% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(hObject);
%
% --------------------------  end  ---------------------------------------