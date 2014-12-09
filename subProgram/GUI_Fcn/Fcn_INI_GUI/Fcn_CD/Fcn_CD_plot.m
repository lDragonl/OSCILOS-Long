function Fcn_CD_plot(hAxes,handles,indexLegend)
% This function is used to plot the schematic view of the combustor
% The input are the tag of axes and handles of current gui window
% first created: 2014-12-03
%
global CI
hFontsize1  = handles.FontSize(1);
poshAxes    = get(handles.axes1,'position');
% -------------------------------------------------------------------------
%
x_sample =      CI.CD.x_sample;                
r_sample =      CI.CD.r_sample;
index    =      CI.CD.index;
%-------------------------------------
W           = abs(x_sample(end) - x_sample(1));             % Length of the combustor
H           = 2*max(r_sample);                              % Diameter of the combustor
plot_ratio  = 1.1;                                          % Ratio of the axes limit to the combustor dimension
axes_W      = plot_ratio*W;                                 % axes x width
axes_H      = 2*plot_ratio*H;                               % axes y width        
x_min       = x_sample(1) - (axes_W-W)./2;                  % axes x min
y_min       = -axes_H./2;                                   % axes y min
%--------------------------------------
cla(hAxes);                             % clear the current axes 
axis(hAxes);
hold on
% -------------------------------------
% These are used for the legend
xMax = max(get(hAxes,'xlim'));
yMax = max(get(hAxes,'ylim'));
plot(hAxes,xMax,yMax,'-b','linewidth',3)
plot(hAxes,xMax,yMax,'-g','linewidth',3)
plot(hAxes,xMax,yMax,'-r','linewidth',3)
plot(hAxes,xMax,yMax,'-m','linewidth',3)
%--------------------------------------
% plot the approximate profile of the combustor which consisting of several
% sections 
for s = 1:length(x_sample)-1
    plot(hAxes,1e3*[x_sample(s),x_sample(s+1)],1e3*[r_sample(s),r_sample(s)],'-k','linewidth',2);
    plot(hAxes,1e3*[x_sample(s),x_sample(s+1)],-1e3*[r_sample(s),r_sample(s)],'-k','linewidth',2);
end
for s = 1:length(x_sample)-2
    plot(hAxes,1e3*[x_sample(s+1),x_sample(s+1)],1e3*[r_sample(s),r_sample(s+1)],'-k','linewidth',2);
    plot(hAxes,1e3*[x_sample(s+1),x_sample(s+1)],-1e3*[r_sample(s),r_sample(s+1)],'-k','linewidth',2);
end
% ---
for s = 1:length(index)
    switch index(s)     % define the color of different sectional interfaces
        case 0          % just interface
            indexColor = 'none';
        case 10         % with heat addition but perturbation
            indexColor = 'm';
        case 11         % with heat addition and heat perturbations
            indexColor = 'r';
    end 
    if s == 1
        indexColor = 'b';   % inlet
    elseif s == length(index)
        indexColor = 'g';   % outlet
    end
    % use a narrow rectangle to represent the type of the interface
    rectangle('parent',hAxes,'Position',1000*[x_sample(s)-axes_W/400,-r_sample(s),...
                axes_W/200,...
                2*r_sample(s)],...
                'Curvature',[0,0],'LineWidth',1,'LineStyle','-',...
                'facecolor',indexColor,'edgecolor',indexColor);
end
set(hAxes,'xlim',1000*[x_min, x_min+axes_W]);
set(hAxes,'box','on','linewidth',0.5,'gridlinestyle','-.');
set(hAxes,'xgrid','on','ygrid','on');
set(hAxes,'ylim',1000*[y_min, y_min+axes_H]);
xlabel(hAxes,'x~ [mm]','Color','k','Interpreter','LaTex');
ylabel(hAxes,'r~ [mm]','Color','k','Interpreter','LaTex');   
%
switch indexLegend
    case 1
    newline = char(10);
    legend1 = ['with mean heat addition', newline, 'and heat perturbations'];
    legend2 = ['with mean heat addition', newline, 'but heat perturbations'];
    hlegend = legend(hAxes,...
                            'inlet','outlet',...
                            legend1,...
                            legend2);
    set(hlegend,'fontsize',hFontsize1,'location','northeastoutside');
    otherwise
end
set(handles.axes1,      'units', 'points',...
                        'Fontunits','points',...
                        'position',poshAxes); 
%
% ------------------------------end----------------------------------------              