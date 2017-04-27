clearvars
%close all

fact = 2;

L = 1;
Lf = 0.1*L;
r = 0.017667350476796;
%r = 0.05;
Tratio = 2;
Mach_num = 0.01*fact;
Flame_Speed = 0.4027*fact;
%Ugs = 4.576;
%Convective_Speed = 5.7206;
fMin = 0;
fMax = 500;
Nfreq = 40;
sigThr = 100;
Nsig = 9;

global CI

for Convective_Speed = 5.7206*fact;
    
    disp(num2str(Convective_Speed))
    CI = [];
    
    %% Set up constants
    CI.Ru       = 8.3145;
    CI.W_air    = 28.96512;
    CI.R_air    = CI.Ru./CI.W_air.*1000;
    
    CI.SD.name_program='subProgram';          % sub-program routine
    addpath(genpath('./'))                    % add directories to search path
    
    % Initialization pannels
    CI.IsRun.GUI_INI_CD         = 0;
    CI.IsRun.GUI_INI_PD         = 0;
    CI.IsRun.GUI_INI_TP         = 0;
    CI.IsRun.GUI_INI_FM_Sel     = 0;
    CI.IsRun.GUI_INI_FM         = 0;
    CI.IsRun.GUI_INI_GEQU       = 0;
    CI.IsRun.GUI_INI_FMEXP      = 0;
    CI.IsRun.GUI_INI_BC         = 0;
    CI.IsRun.GUI_INI_BC_Entropy = 0;
    
    % Frequency domain calculation pannels
    CI.IsRun.GUI_FREQ_EigCal    = 0;
    CI.IsRun.GUI_FREQ_EigCal_AD = 0;
    
    % Time domain simulation
    CI.IsRun.GUI_TD_Convg       = 0;
    CI.IsRun.GUI_TD_Para_Config = 0;
    CI.IsRun.GUI_TD_Cal_JLI_AMorgans = 0;
    CI.IsRun.GUI_TD_Cal_simple = 0;
    
    %% Set up Geometry
    CI.CD.x_sample = [0 Lf L];
    CI.CD.r_sample = [r r r];
    CI.CD.SectionIndex = [0 11 0];
    CI.CD.TubeIndex = [0 0 0];
    CI.CD.pop_CD_type = 1;
    CI.CD.NUM_HR = 0;
    CI.CD.HR_config = [];
    CI.CD.NUM_Liner = 0;
    CI.CD.Liner_config = [];
    CI.CD.indexHA = 2;
    CI.CD.indexHP = 2;
    CI.CD.indexLiner =  find(CI.CD.SectionIndex == 30);
    CI.CD.indexHR = find(CI.CD.SectionIndex == 2);
    CI.CD.isHA = 1;
    CI.CD.isHP = 1;
    CI.CD.isLiner = 0;
    CI.CD.isHR = 0;
    
    %% Set up Thermal Properties
    gamma_cst = 1.4;
    Cp_cst = 1005;
    
    CI.TP.indexHA = 2;
    CI.TP.HA_style = 1;
    CI.TP.TRatio = Tratio;
    CI.TP.indexFuel = 1;
    CI.TP.eff = 1;
    CI.TP.Phi = 1;
    CI.TP.Q = 0;
    CI.TP.DeltaHr = 0;
    CI.TP.HA_num = 1;
    
    CI.TP.numSection        = length(CI.CD.x_sample)-1;
    CI.TP.p_mean(1:2,1)     = 101325;
    CI.TP.T_mean(1:2,1)     = 293.15;
    
    CI.TP.index_gamma = 2;
    CI.TP.c_mean(1:2,1)     = (gamma_cst.*CI.R_air.*CI.TP.T_mean(1:2,1)).^0.5;
    CI.TP.Cp(1:2,1)         = Cp_cst;
    CI.TP.gamma(1:2,1)      = gamma_cst;
    
    CI.TP.M_mean(1:2,1)     = Mach_num;
    CI.TP.u_mean(1:2,1)     = CI.TP.M_mean(1,1).*CI.TP.c_mean(1,1);      % mean velocity
    CI.TP.rho_mean(1:2,1)   = CI.TP.p_mean(1,1)./(CI.R_air.*CI.TP.T_mean(1,1)); % mean density
    
    for ss = 1:CI.TP.numSection-1
        CI.TP.Theta(ss)  = (CI.CD.r_sample(ss+1)./CI.CD.r_sample(ss)).^2;       % Surface area ratio S2/S1
    end
    
    % Begin adding by Dong Yang
    Liner_Flag = 1; % Initialize the flag of the Liner considered by the loop
    HR_Flag = 1;    % Initialize the flag of the HR considered by the loop
    % End adding by Dong Yang
    % --------------------------------
    indexHA_num = 0;                % set the initial value to zero and it will be increased by 1 after every HA interface
    for ss = 1:CI.TP.numSection-1
        % In every interface, the changes are splitted to two steps:
        % 1. cross sectional surface area change
        % 2. HA or .....
        % --------------step 1-------------------------
        %
        [   CI.TP.p_mean(1:2,ss+1),...
            CI.TP.rho_mean(1:2,ss+1),...
            CI.TP.u_mean(1:2,ss+1) ]...
            = Fcn_calculation_TP_mean_WO_HeatAddition(CI.TP.p_mean(1,ss),...
            CI.TP.rho_mean(1,ss),...
            CI.TP.u_mean(1,ss),...
            CI.TP.Theta(ss),...
            CI.TP.gamma(1,ss));
        % ----------
        CI.TP.gamma(1:2,ss+1)   = CI.TP.gamma(1,ss);
        CI.TP.Cp(1:2,ss+1)      = CI.TP.Cp(1,ss);
        CI.TP.T_mean(1:2,ss+1)  = CI.TP.gamma(1,ss+1)/(CI.TP.gamma(1,ss+1)-1)...
            *CI.TP.p_mean(1,ss+1)./(CI.TP.Cp(1,ss+1).*CI.TP.rho_mean(1,ss+1));
        CI.TP.c_mean(1:2,ss+1)  = ((CI.TP.gamma(1,ss+1) - 1).*CI.TP.Cp(1,ss+1).*CI.TP.T_mean(1,ss+1)).^0.5;
        CI.TP.M_mean(1:2,ss+1)  = CI.TP.u_mean(1,ss+1)./CI.TP.c_mean(1,ss+1);
        % ----------
        %
        % --------------step 2-------------------------
        %
        switch CI.CD.SectionIndex(ss+1)
            case 0
                % in case 0, no changes
            case {10,11}                                    % with HA
                indexHA_num = indexHA_num + 1;              % this number is increased by 1
                % ---------first calculate the temperature, Cp, Rg after the
                % heat addition interface -------------------------------------
                [   CI.TP.TRatio(indexHA_num),...
                    CI.TP.c_mean(1,ss+1),...
                    CI.TP.DeltaHr(indexHA_num),...          % heat release rate per mass flow rate
                    CI.TP.gamma(1,ss+1),...
                    CI.TP.Cp(1,ss+1)] =...
                    Fcn_GUI_INI_TP_calculation_products_after_HA(   [],...
                    indexHA_num,...
                    CI.TP.T_mean(2,ss+1),...
                    CI.TP.p_mean(2,ss+1));
                CI.TP.T_mean(1,ss+1)    = CI.TP.TRatio(indexHA_num).*CI.TP.T_mean(2,ss+1);
                switch CI.TP.index_gamma
                    case 1
                        CI.TP.gamma(1,ss+1)     = gamma_cst;
                        CI.TP.Cp(1,ss+1)        = Cp_cst;
                        CI.TP.c_mean(1,ss+1)    = (gamma_cst*CI.R_air*CI.TP.T_mean(1,ss+1)).^0.5;
                    case 2
                        % nothing happens!
                end
                Rg2 = CI.TP.Cp(1,ss+1)./(CI.TP.gamma(1,ss+1)./(CI.TP.gamma(1,ss+1)-1));
                % ---------then, use the resolved temperature, Rg and the mean
                % properties after the area changes to calculate the mean
                % properties after HA ----------------------------------------
                [   CI.TP.p_mean(1,ss+1),...
                    CI.TP.rho_mean(1,ss+1),...
                    CI.TP.u_mean(1,ss+1)] = ...
                    Fcn_calculation_TP_mean_W_HeatAddition( CI.TP.p_mean(2,ss+1),...
                    CI.TP.rho_mean(2,ss+1),...
                    CI.TP.u_mean(2,ss+1),...
                    Rg2,...
                    CI.TP.T_mean(1,ss+1));
                CI.TP.M_mean(1,ss+1)  = CI.TP.u_mean(1,ss+1)./CI.TP.c_mean(1,ss+1);
                switch CI.TP.index_gamma
                    case 1  % constant Cp and gamma
                        CI.TP.DeltaHr(indexHA_num) = CI.TP.Cp(1,ss+1).*(CI.TP.T_mean(1,ss+1) - CI.TP.T_mean(2,ss+1))...
                            + 0.5*(CI.TP.u_mean(1,ss+1).^2 - CI.TP.u_mean(1,ss).^2);
                    case 2
                        % nothing happens!
                end
                mass = CI.TP.rho_mean(2,ss+1).*CI.TP.u_mean(2,ss+1).*CI.CD.r_sample(ss+1).^2.*pi; % mass flow rate before HA
                CI.TP.Q(indexHA_num)  = CI.TP.DeltaHr(indexHA_num).*mass;       % heat release rate
                
                % Begin adding by Dong Yang
            case 2
                % Typically, it can be assumed that mean flow will not be affect by the HR, anyway,
                % the following procedure is an alternate way to calculate mean parameters across the cross section where the resonator is installed.
                
                section_Num=ss;
                [   CI.TP.T_mean(1,ss+1),...
                    CI.TP.rho_mean(1,ss+1),...
                    CI.TP.u_mean(1,ss+1),...
                    CI.TP.Cp(1,ss+1),...
                    CI.TP.gamma(1,ss+1)] = ...
                    Fcn_calculation_TP_mean_across_HR( HR_Flag,...
                    section_Num,...
                    CI.TP.T_mean(2,ss+1),...
                    CI.TP.rho_mean(2,ss+1),...
                    CI.TP.u_mean(2,ss+1));
                CI.TP.p_mean(1:2,ss+1)  = CI.TP.T_mean(1,ss+1).*CI.TP.rho_mean(1,ss+1)*(CI.TP.gamma(1,ss+1)-1)./CI.TP.gamma(1,ss+1)*CI.TP.Cp(1,ss+1);
                CI.TP.c_mean(1:2,ss+1)  = ((CI.TP.gamma(1,ss+1) - 1).*CI.TP.Cp(1,ss+1).*CI.TP.T_mean(1,ss+1)).^0.5;
                CI.TP.M_mean(1:2,ss+1)  = CI.TP.u_mean(1,ss+1)./CI.TP.c_mean(1,ss+1);
                HR_Flag = HR_Flag+1;
            case 30
                % Nothing happens because mean flow properties are constant
                % across the inlet side cross-section of the lined duct
            case 31
                % It is assumed only mean flow velocity will increase a bit in
                % the lined duct region.
                [ CI.TP.u_mean(1,ss+1)] = Fcn_calculation_TP_mean_across_Liner(Liner_Flag,ss+1);
                CI.TP.M_mean(1,ss+1)=CI.TP.u_mean(1,ss+1)/CI.TP.c_mean(1,ss+1);
                Liner_Flag = Liner_Flag+1;
                % End adding by Dong Yang
        end
    end
    
    
    %% Set up Flame model
    NumHP = length(CI.CD.indexHP);
    CI.FM.ModelType = {     'Linear FTF model';...
        'Nonlinear FDF model';...
        'Experimental/CFD fitted FDF';...
        'Fully non-linear G-Equation model';...
        'Fully nonlinear convective G-equation model'};
    for ss = 1:NumHP
        CI.FM.indexFM(ss)       =   5;      % defult setting: linear FTF
        CI.FM.HP{ss}         =   [];
    end   
    
    CI.FM.HP{NumHP}.NL.style      = 5; % G-Equation (Markstein 1964)
    CI.FM.HP{NumHP}.GEQU_CONV.nb_points = 400; % Number of points used for discretisation along r
    CI.FM.HP{NumHP}.GEQU_CONV.nf_points = 400; % Number of points used for discretisation along r
    % AO: To me it is a bad idea to treat this as a vecotr. It forces all the
    % flame moels to be the same. Just set every flame indipendently...
    CI.FM.HP{NumHP}.GEQU_CONV.rb = CI.CD.r_sample(CI.CD.indexHP); % If there are multiple heat zones in the duct, this is a vector
    CI.FM.HP{NumHP}.GEQU_CONV.ra   = CI.FM.HP{NumHP}.GEQU_CONV.rb/2; % in m, also a vector
    CI.FM.HP{NumHP}.GEQU_CONV.area_ratio = 1.0 -(CI.FM.HP{NumHP}.GEQU_CONV.ra./CI.FM.HP{NumHP}.GEQU_CONV.rb).^2;
    CI.FM.HP{NumHP}.GEQU_CONV.rho1 = CI.TP.rho_mean(1,max(CI.CD.indexHP - 1,1));% This is a vector if there are multple flame. The max function is required is the flame is at the begining of the duct.
    CI.FM.HP{NumHP}.GEQU_CONV.U1 = CI.TP.u_mean(1,max(CI.CD.indexHP - 1,1));% This is a vector if there are multple flame. The max function is required is the flame is at the begining of the duct.
    CI.FM.HP{NumHP}.GEQU_CONV.Ugs   = CI.FM.HP{NumHP}.GEQU_CONV.U1/CI.FM.HP{NumHP}.GEQU_CONV.area_ratio; % in m/s, this is a vector if there are multple flame.
    CI.FM.HP{NumHP}.GEQU_CONV.SU   = Flame_Speed; % in m/s, this is a vector if there are multple flame.
    CI.FM.HP{NumHP}.GEQU_CONV.UC   = Convective_Speed; % in m/s, this is a vector if there are multple flame.
    CI.FM.HP{NumHP}.GEQU_CONV.CFL   = 0.5; % in m/s, this is a vector if there are multple flame.
    CI.FM.HP{NumHP}.GEQU_CONV.F_Cutoff = 750; % in m/s, this is a vector if there are multple flame.
    
    for runner = 1:length(CI.FM.HP{NumHP}.GEQU_CONV.ra)
        CI.FM.HP{NumHP}.GEQU_CONV.y_vec(runner,:) = linspace(CI.FM.HP{NumHP}.GEQU_CONV.ra(runner),CI.FM.HP{NumHP}.GEQU_CONV.rb(runner),...
            CI.FM.HP{NumHP}.GEQU_CONV.nb_points(runner)); % currently the nb of points for all flames nees to be the same
    end
    
    CI.FM.HP{NumHP}.GEQU_CONV.Lf = sqrt((CI.FM.HP{NumHP}.GEQU_CONV.Ugs/CI.FM.HP{NumHP}.GEQU_CONV.SU)^2 - 1)*...
        (CI.FM.HP{NumHP}.GEQU_CONV.rb-CI.FM.HP{NumHP}.GEQU_CONV.ra); % Flame length
    CI.FM.HP{NumHP}.GEQU_CONV.dr = 2*CI.FM.HP{NumHP}.GEQU_CONV.Lf/CI.FM.HP{NumHP}.GEQU_CONV.nf_points;  %dr axial
    [CI.FM.HP{NumHP}.GEQU_CONV.xi_steady, CI.FM.HP{NumHP}.GEQU_CONV.x_flow] = ...
        Fcn_TD_Conv_Gequ_steady_flame( CI.FM.HP{NumHP}.GEQU_CONV.Ugs,CI.FM.HP{NumHP}.GEQU_CONV.SU,CI.FM.HP{NumHP}.GEQU_CONV.y_vec,...
        CI.FM.HP{NumHP}.GEQU_CONV.dr );
    CI.FM.HP{NumHP}.GEQU_CONV.uflow = zeros(size(CI.FM.HP{NumHP}.GEQU_CONV.x_flow));
    CI.FM.HP{NumHP}.GEQU_CONV.xi = CI.FM.HP{NumHP}.GEQU_CONV.xi_steady; % initialise xi value
    
    % Set the flame transfer function numerator and denominator to 1, as they are used for the Green's function initilisation
    CI.FM.HP{NumHP}.FTF.num = 1;
    CI.FM.HP{NumHP}.FTF.den = 1;
    
    
    %% Boundary Conditions
    
    CI.BC.StyleInlet = 2;    
    CI.BC.StyleOutlet = 1;
    
    CI.BC.num1 = 1;
    CI.BC.den1 = 1;
    CI.BC.tau_d1 = 0;
    
    CI.BC.num2 = -1;
    CI.BC.den2 = 1;
    CI.BC.tau_d2 = 0;
    
    CI.BC.A1 = 1;
    CI.BC.Phi1 = 0;
    
    CI.BC.A2 = 1;
    CI.BC.Phi2 = 0;
    
    CI.BC.ET = [];
    CI.BC.ET.pop_type_model = 1;
    CI.BC.ET.Dispersion.Delta_tauCs = 0;
    CI.BC.ET.Dissipation.k=0;
    
    %% Calculate EigenFrequencies
    %function GUI_FREQ_EigCal('OSCILOS_long', handle_figure);
    Fcn_PreProcessing
    
    CI.EIG.APP_style = 13;
    CI.EIG.Scan = [];
    CI.EIG.Scan.FreqMin = fMin;
    CI.EIG.Scan.FreqMax = fMax;
    CI.EIG.Scan.FreqNum = Nfreq;
    CI.EIG.Scan.GRMin   = -sigThr;
    CI.EIG.Scan.GRMax   = sigThr;
    CI.EIG.Scan.GRNum   = Nsig;
    
    for ss = 1:length(CI.CD.indexHP)
        CI.FM.indexHPinHA(ss) = find(CI.CD.indexHA == CI.CD.indexHP(ss));
    end
    
    %GUI_FREQ_EigCal_Eigenvalues_calculation(hObject);
    GUI_FREQ_EigCal_Scan_domain_Contour_domain_splitting;
    
    %function Fcn_Para_initialization
    CI.EIG.FDF.uRatioNum    = 1;
    CI.EIG.FDF.uRatioSp(1)  = 0;
    CI.EIG.FDF.pRatioNum    = 1;
    CI.EIG.FDF.pRatioSp(1)  = 0;
    CI.EIG.FDF.num{1}       = [];
    CI.EIG.FDF.den{1}       = [];
    CI.EIG.FDF.tauf(1)      = 0;
    
    % Fcn_calculation_APP_1(hObject)
    CI.EIG.Scan.EigValCol{ss}   = Fcn_calculation_eigenvalues(1);
    CI.EIG.Cont.ValCol{ss}      = Fcn_calculation_contour(1);
    
    eigenvalue = CI.EIG.Scan.EigValCol{1};
    
    %% Plotting
    indexShow = 1;
    CI.EIG.pop_numMode = 1;
    Eigenvalue          = CI.EIG.Scan.EigValCol{indexShow};
    ValueContour        = CI.EIG.Cont.ValCol{indexShow};
    
    figure(1)
    clf
    hAxes1 = gca;
    hAxes2 = gca;
    fontSize1 = 10;
    fontSize2 = 10;
    
    hold on
    %
    contourf(hAxes1,CI.EIG.Cont.GRSp./100,CI.EIG.Cont.FreqSp,20*log10(abs(ValueContour')))
    drawnow
    position_hAxes1=get(hAxes1,'position');
    ylimitUD = [CI.EIG.Scan.FreqMin CI.EIG.Scan.FreqMax];
    xlimitUD = [CI.EIG.Scan.GRMin   CI.EIG.Scan.GRMax]./100;
    hold off
    %
    set(hAxes1,'YColor','k','Box','on');
    set(hAxes1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
    xlabel(hAxes1,  '$ Re(s)/100: \textrm{Growth rate}~~/100~~$ [rad s$^{-1}$] ',...
        'Color','k','Interpreter','LaTex','FontSize',fontSize1);
    ylabel(hAxes1,'$ Im(s)/2\pi: \textrm{Frequency}~~$ [Hz]','Color','k',...
        'Interpreter','LaTex','FontSize',fontSize1);
    set(hAxes1,'ylim',ylimitUD,'YAxisLocation','left','Color','w');
    set(hAxes1,'xlim',xlimitUD);
    grid on
    colorbar
    colormap(hot);
    hcb=colorbar;
    set(hcb,'Fontsize',fontSize2,'box','on','Unit','points')
    set(hcb,'position',[position_hAxes1(1)+position_hAxes1(3),...
        position_hAxes1(2),...
        position_hAxes1(3)./20,...
        position_hAxes1(4).*1]);
    set(hAxes1,'position',position_hAxes1)
    hcb_ylim=get(hcb,'ylim');
    handles.hColorbar.ylimit=[  min(min(20*log10(abs(ValueContour')))),...
        max(max(20*log10(abs(ValueContour'))))];
    %------------------------------------
    
    axes(hAxes2)
    set(hAxes2,'position',get(hAxes1,'position'));
    hold on
    plot(hAxes2,real(Eigenvalue)./100,imag(Eigenvalue)./2./pi,'p',...
        'markersize',8,'color','k','markerfacecolor',[1,1,1])
    drawnow
    hold off
    set(hAxes2,     'ylim', get(hAxes1,'ylim'),...
        'yTick',get(hAxes1,'ytick'),...
        'YAxisLocation','left','Color','none');
    set(hAxes2,     'xlim', get(hAxes1,'xlim'),...
        'xTick',get(hAxes1,'xtick'),...
        'xcolor','b','ycolor','b','gridlinestyle','-.');
    set(hAxes2,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',0.5)
    grid on
    %------------------------
    xt_pos=min((get(hAxes2,'xlim')))+1.15*(max((get(hAxes1,'xlim'))) - min((get(hAxes1,'xlim'))));
    yt_pos=mean(get(hAxes2,'ylim'));
    hTitle = title(hAxes2, 'Eigenvalues are located at minima');
    set(hTitle, 'interpreter','latex', 'fontunits','points','fontsize',fontSize1)
    
    hold on
    plot([0 0],[0 CI.EIG.Scan.FreqMax],'k--','linewidth',2)
    hold off
    drawnow
    % % Eigenmodes
    %     s_star = Eigenvalue(CI.EIG.pop_numMode);             % eigenvalue
    %         [x_resample,p,u] = Fcn_calculation_eigenmode_Linear(s_star);
    
end

