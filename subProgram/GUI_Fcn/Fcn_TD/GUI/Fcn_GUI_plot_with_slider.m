function f = Fcn_GUI_plot_with_slider(data,slider_vec,x_vec,titles,fig,axes)

% This function creates a plot of a matrix with a slider. 
% the data matrix should have the x_values as columns, and the slider
% values as lines

pannelsize=get(axes,'position');
pW=pannelsize(3);
pH=pannelsize(4);
% set(handles.axes1,      'units', 'points',...
%     'Fontunits','points',...
%     'position',[pW*2/10 pH*2.0/10 pW*7.5/10 pH*6.5/10],...
%     'fontsize',handles.FontSize(1),...
%     'color',handles.bgcolor{1},...
%     'box','on');

% Create plot
f = fig;
h = plot(x_vec,data(1,:));

if (nargin == 4) % labels have been provided for the functions
    xlabel(titles(1))
    ylabel(titles(2))
end

% add the slider
b = uicontrol('Parent',f,'Style','slider','Position',[pW*2/10 pH*0.5/10 pW*7.5/10 pH*0.5/10],...
     'value',1, 'min',1, 'max',length(slider_vec));
bgcolor = get(f,'Color'); % colour of the figure
bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
    'String',num2str(slider_vec(1)),'BackgroundColor',bgcolor);
bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
    'String',num2str(slider_vec(end)),'BackgroundColor',bgcolor);
bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
    'String',[titles(3),' = ',num2str(slider_vec(1))],'BackgroundColor',bgcolor);

% add a listener for the slider handle b
try    % R2013b and older, if it is still valid and you get a warning to use ContinuousValueChange, disregard warning
    addlistener(b,'ActionEvent',@(src,eventdata)slider_callback(src,eventdata,data,slider_vec,h,f,bgcolor,titles(3)));
catch  % R2014a and newer
    addlistener(b,'ContinuousValueChange',@(src,eventdata)slider_callback(src,eventdata,data,slider_vec,h,f,bgcolor,titles(3)));
end

end

function slider_callback(src,evt,data,slider_vec,plot_handle,figure_handle,bgcolor,slider_title)
slider_value=round(get(src,'value'));
set(plot_handle,'YData',data(slider_value,:))
uicontrol('Parent',figure_handle,'Style','text','Position',[240,25,100,23],...
    'String',[titles(3),' = ',num2str(slider_vec(slider_value))],'BackgroundColor',bgcolor);
end