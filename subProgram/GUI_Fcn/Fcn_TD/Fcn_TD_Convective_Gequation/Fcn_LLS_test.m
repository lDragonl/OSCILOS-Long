function [] = Fcn_LLS_test(A)
%
addpath(genpath('./'))           % Add directories to search path
%
%   A = 0.05;
    parfor fps = 1:1:20
        fp          = fps*20;               % Frequency of perturbations
        %
        flameSpeed  = 0.368;             % Flame Speed
        vMean       = 1.5;               % Mean Flow Speed
        vConv       = 1.5;               % Convective Speed
        ra          = 25/1000;           % Internal radius
        rb          = 50/1000;           % External radius
        Lf          = sqrt(vMean^2/0.368^2-1)*(rb-ra);
        
        dx1         = (rb-ra)/50;
        dx2         = vMean/fp/30;
        dx          = min(dx1,dx2);
        dx          = ra/(round(ra/dx));
        dy          = dx;
        
        CFL         = 0.3;               % CFL number
        dt          = dx*CFL/vConv;      % Time step
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [dataGhost, gridGhost, SET, tNow, vyPrime_no_stencil] = FcnLLS_Main_INI(vMean, vConv, flameSpeed, dx, dy, dt, ra, rb);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        gridGhost.RA = ra;
        %NumStep     = 2000;
        NumStep = round((Lf/vMean + 3/fp)/dt)
        %
        NGapPlot    = round(1/fp/100/SET.dt);
        %
        tSp         = tNow : SET.dt : (NumStep)*SET.dt;
        
        tEnd        = (NumStep)*SET.dt;
        % chirp signal
        % vyPrime_inlet_Sp = 0.8*sin(2*pi*(fp + fp.*tSp./tEnd).*tSp);
        % sine signal
        vyPrime_inlet_Sp = A*sin(2*pi*(fp).*tSp);
        %
        idOld       = 1;
        QRatio      = zeros(1, NumStep);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % figure
        hff = figure;
        scrsz = get(0,'ScreenSize');
        set(hff,'Position',[scrsz(4).*(1/8) scrsz(4).*(1/20) scrsz(3)*0.75 scrsz(4)/3*2])
        hAxes1=axes('Unit','pixels','position',[100 100 300 450]);
        
        hAxes2=axes('Unit','pixels','position',[500 100 500 450]);
        %
        pause(1)
        for ss = (1 : NumStep)
            Niter = ss;
            
            %        tic
            if(mod(ss,100)==0)
                disp([num2str(ss/NumStep*100,'%.1f'),'%']);
            end
            isReINI_Global = 0;
            if rem(Niter, 10) == 0
                isReINI_Global = 1;
            end
            [dataGhost, FF, vyPrime_no_stencil, tNow] = ...
                FcnLLS_MainCal_Dynamic(tNow, SET, dataGhost, gridGhost, vyPrime_no_stencil, vyPrime_inlet_Sp(ss + 1), isReINI_Global);
            %
            if ss == 1
                FF0 = FF;
            end
            QRatio(ss) = FF.A./FF0.A;
            %
            %     idOld       = Fcn_FF_Plot(hAxes1, hAxes2, FF, idOld, ss, NGapPlot, FF0);
            %         if(mod(ss,20)==0)
            %             figure(1)
            %             clf
            %             contourf(gridGhost.xs{1},gridGhost.xs{2},dataGhost,100,'LineColor','none')
            %             hold on
            %             contour(gridGhost.xs{1},gridGhost.xs{2},dataGhost,[0 0],'LineColor','k')
            %             for ii = 1:length(FF.xFFInterp)
            %                 plot(FF.xFFInterp{ii}, FF.yFFInterp{ii}, '.', 'markersize', 8, 'color', 'r')
            %             end
            %             caxis([min(min(dataGhost)) max(max(dataGhost))])
            %             pause(0.1)
            %             drawnow
            %         end
            
            %       toc
        end
        % save('chirp_4times.mat', 'tSp', 'vyPrime_inlet_Sp', 'QRatio');
        
        Sfp{fps} = fp;
        StSP{fps} = tSp;
        SvyPrime_inlet_Sp{fps} = vyPrime_inlet_Sp;
        SQRatio{fps} = QRatio;
    end
    
    
    
    for ii=1:length(Sfp)
        filename         = ['FDF/FDF_A',num2str(A,'%.2f'),'_f_',num2str(Sfp{ii}),'.mat'];
        tSp              = StSP{ii};
        vyPrime_inlet_Sp = SvyPrime_inlet_Sp{ii};
        QRatio           = SQRatio{ii};
        save(filename, 'tSp', 'vyPrime_inlet_Sp', 'QRatio');
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  breaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function idOld = Fcn_FF_Plot(hAxes1, hAxes2, FF, idOld, idNew, NGapPlot, FF0)
if idNew == idOld + NGapPlot
    hold(hAxes1, 'off')
    for ii = 1:length(FF.xFFInterp)
        plot(hAxes1, FF.xFFInterp{ii}, FF.yFFInterp{ii}, '.', 'markersize', 8, 'color', 'r')
        hold(hAxes1, 'on')
        plot(hAxes1, -FF.xFFInterp{ii}, FF.yFFInterp{ii}, '.', 'markersize', 8, 'color', 'r')
    end
    hold(hAxes1, 'off')
    ylim(hAxes1,[0, 0.12])
    xlabel(hAxes1,'x [m]', 'Fontsize', 15)
    ylabel(hAxes1,'y [m]', 'Fontsize', 15)
    grid(hAxes1, 'on')
    title(hAxes1,'Flame front evolution', 'Fontsize', 15)
    set(hAxes1, 'Fontsize', 15)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    hold(hAxes2, 'on')
    plot(hAxes2, 1e3.*FF.tNow, FF.A./FF0.A, '.', 'markersize', 8, 'color', 'r')
    %     hold on
    %     plot(hAxes2, 1e3.*FF.tNow, FF.Anew./FF0.Anew, '.', 'markersize', 8, 'color', 'k')
    xlabel(hAxes2,'t [ms]','Fontsize', 15)
    ylabel(hAxes2,'Qratio','Fontsize', 15)
    grid(hAxes2, 'on')
    title(hAxes2,'Evolution of normalised heat release rate perturbation','Fontsize', 15)
    set(hAxes2, 'Fontsize', 15)
    %%%%
    pause(0.1)
    idOld = idNew;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  breaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



