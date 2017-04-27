function CI = Fcn_Eig_WaveBased
% last updated on 2017-02-09, by Jingxuan Li
%
close all
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(genpath('./'))                    % add directories to search path
% rmpath(genpath('./'))                    % remove directories from search path
%%%%%%%%%%%%%%%%%%%%%%%%%%% Main part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

Tjump = linspace(1.01,10,100);
for ii = 1:length(Tjump)
    CI          = Fcn_INI('Tjump',Tjump(ii));
    %
    CaseName    = ['Case_1'];
    Fcn_ContourMap(@MyFunDet_NoFlow, CI, CaseName);
    %
    assignin('base','CI',CI);
    
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Main part %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function CI = Fcn_INI(parType,par)
%
CI = Fcn_CD; %Geometry
CI = Fcn_TP(CI,parType,par);
CI = Fcn_FM(CI);
CI = Fcn_BC(CI);
CI = Fcn_PreProcessing(CI);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Fcn_ContourMap(myFun, CI, CaseName)
%
FreqSpCM    = 20:10:600;
GrSpCM      = -800:10:800;
%
F           = @(s, CI)myFun(s, CI);
Val         = Fcn_ValCal_Contour_Map(F, CI, FreqSpCM, GrSpCM);
%
Freqfsolve  = linspace(FreqSpCM(1),    FreqSpCM(end), 10);
Grfsolve    = linspace(GrSpCM(1),      GrSpCM(end),   10);
%
EigVal      = Fcn_EigCal(F, CI, Freqfsolve, Grfsolve);
%
Fcn_ContourMap_Plot(Val, EigVal, FreqSpCM, GrSpCM, CaseName);
%
EIG.FreqSpCM    = FreqSpCM;
EIG.GrSpCM      = GrSpCM;
EIG.Freqfsolve  = Freqfsolve;
EIG.Grfsolve    = Grfsolve;
EIG.EigVal      = EigVal;
EIG.Val         = Val;
%
MatName = ['Mat_ContourMap_' CaseName '.mat'];
save(MatName, 'EIG');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Val = Fcn_ValCal_Contour_Map(myFun, CI, FreqSp, GrSp)
n1 = length(FreqSp);
n2 = length(GrSp);
Val = zeros(n1, n2);
for ss = 1 : n1
    for kk = 1 : n2
        s = GrSp(kk) + 1i*2*pi*FreqSp(ss);
        F   = @(s, CI)myFun(s, CI);
        Val(ss, kk) = feval(F, s, CI);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function EigVal = Fcn_EigCal(myFun, CI, FreqSp, GrSp)
F = @(s, CI)myFun(s, CI);
EigVal    = FcnCalEigVal1(F,       FreqSp, GrSp, CI);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function Fcn_ContourMap_Plot(Val, EigVal, FreqSpCM, GrSpCM, CaseName)
%%%%%%%%%%%%%%%%%%
% FigName = 'Figure';
% idExist = isdir(FigName);
% if idExist==0
%     mkdir(FigName);
% end
%%%%%%%%%%%%%%%%%%
isPlot_fsolve = 1;
% --------------------------  PLOT  ---------------------------------------
h=figure
fontSize1 = 20;
fontSize2 = 20;
scrsz = get(0,'ScreenSize');
set(h,'Position',[scrsz(4).*(1/8) scrsz(4).*(1/20) scrsz(3)*2/5 scrsz(4).*(4/5)])
%************
hAxes1_figure1=axes('Unit','pixels','position',[100 100 400 400]);
position_hAxes1_figure1=get(hAxes1_figure1,'position');
hold on
contourf(GrSpCM./100, FreqSpCM, log10(abs(Val)))
drawnow
hold off
set(hAxes1_figure1,'YColor','k','Box','on');
set(hAxes1_figure1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',1)
h1x_figure1 = xlabel(hAxes1_figure1,  '$ Re(s)/100: \textrm{Growth rate}~~/100~~$ [s$^{-1}$] ',...
    'Color','k','Interpreter','LaTex','FontSize',fontSize2);
h1y_figure1 = ylabel(hAxes1_figure1,'$ Im(s)/2\pi: \textrm{Frequency}~~$ [Hz]','Color','k',...
    'Interpreter','LaTex','FontSize',fontSize2);
set(h1y_figure1,'Unit','pixels','position',[-50 200 0],'visible','on');
set(h1x_figure1,'Unit','pixels','position',[200 -40 0],'visible','on');
grid on
colormap(hot);
hcb_figure1=colorbar;
set(hcb_figure1,'Fontsize',fontSize1,'box','on','Unit','pixels')
set(hcb_figure1,'position',[position_hAxes1_figure1(1),...
    position_hAxes1_figure1(2)+position_hAxes1_figure1(4),...
    position_hAxes1_figure1(3),...
    position_hAxes1_figure1(4)./50],...
    'location','northoutside');
set(hAxes1_figure1,'position',position_hAxes1_figure1)
%************

hAxes2_figure1=axes('Unit','pixels','position',position_hAxes1_figure1);
hold on
switch isPlot_fsolve
    case 1
        plot(hAxes2_figure1,real(EigVal)./100,imag(EigVal)./2./pi,'p',...
            'markersize',10,'color','k','markerfacecolor','w')
        drawnow
        hold off
    otherwise
end
set(hAxes2_figure1,     'ylim', get(hAxes1_figure1,'ylim'),...
    'yTick',get(hAxes1_figure1,'ytick'),...
    'yticklabel',[],...
    'YAxisLocation','left','Color','none');
set(hAxes2_figure1,     'xlim', get(hAxes1_figure1,'xlim'),...
    'xTick',get(hAxes1_figure1,'xtick'),...
    'xticklabel',[],...
    'xcolor','b','ycolor','b','gridlinestyle','-.');
set(hAxes2_figure1,'FontName','Helvetica','FontSize',fontSize1,'LineWidth',0.5)
grid on

%********************
set(gcf,'PaperUnits','points')
set(gcf,'PaperPosition',[1,1,1800,650])
%****************************
savename = ['Fig_' CaseName];
saveas(gcf,savename,'fig');
% saveas(gcf,savename,'epsc')
% savenameEps = [savename '.eps'];
% savenamePdf = [savename '.pdf'];
% eps2pdf(savenameEps,savenamePdf)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CI = Fcn_CD
% coordinates of the combustion dimensions.
% The first index corresponds to the beginning of the combustor
% the final index corresponds to the end of the combustor
rSp = [0.02,   0.02,  0.02];
xSp = [0,     0.3,   1];
%
indexFlame = 2;         % the flame is located at the beginning of module indexFlame
%
N   = length(xSp) - 1;  % total number of modules
%
CI.CD.indexFlame = indexFlame;
CI.CD.xSp   = xSp;
CI.CD.rSp   = rSp;
CI.CD.S     = pi.*rSp.^2;
%
for ss = 1 : N - 1
    CI.CD.Theta(ss) = CI.CD.S(ss + 1) ./ CI.CD.S(ss); % surface area ratio
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
function CI = Fcn_TP(CI,parType,par)
%%%%%%%%%%%%%%%%%
% to calculate the thermal properties in different modules
% the ratio of specific heats is assumed constant
% zero mean flow
%
idf             = CI.CD.indexFlame;
T1              = 300;
T2              = 600;
if(strcmp(parType,'Tjump'))
    T2 = par*T1;
end
p1              = 1e5;
Rg              = 287;
gamma           = 1.4;
CI.TP.Rg        = Rg;
cp              = gamma/(gamma - 1)*Rg;
CI.TP.cp        = cp;
CI.TP.gamma     = gamma;
N               = length(CI.CD.xSp) - 1;
%
CI.TP.pMean(1 : N)  = p1;
CI.TP.TMean(1 : N)  = T2;
CI.TP.TMean(1 : idf - 1)  = T1;
%
CI.TP.cMean         = (gamma .* Rg .* CI.TP.TMean).^0.5;
%
CI.TP.tau(1 : N)    = diff(CI.CD.xSp) ./ (CI.TP.cMean);
CI.TP.rhoMean       = CI.TP.pMean ./ (CI.TP.Rg .* CI.TP.TMean);
%
CI.TP.alpha         = (gamma - 1) .* cp .* (T2 - T1) ./ (CI.TP.cMean(idf - 1).^2);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CI = Fcn_FM(CI)
CI.FM.indexFM = 5;
CI.FM.HP{1}.IsRun = 1;
CI.FM.HP{1}.NL.style = 5;
CI.FM.HP{1}.FTF.num = 1;
CI.FM.HP{1}.FTF.den = 1;
CI.FM.HP{1}.GEQU_CONV.nb_points = 400;
CI.FM.HP{1}.GEQU_CONV.nf_points = 400;
CI.FM.HP{1}.GEQU_CONV.rb = 0.05;
CI.FM.HP{1}.GEQU_CONV.ra = 0.025;
CI.FM.HP{1}.GEQU_CONV.area_ratio = (1-CI.FM.HP{1}.GEQU_CONV.rb^2/CI.FM.HP{1}.GEQU_CONV.ra^2);
CI.FM.HP{1}.GEQU_CONV.rho1 = CI.TP.rhoMean(1);
CI.FM.HP{1}.GEQU_CONV.U1 = CI.TP.uMean(1);
CI.FM.HP{1}.GEQU_CONV.Ugs = CI.TP.rho_mean/CI.FM.HP{1}.GEQU_CONV.area_ratio;
CI.FM.HP{1}.GEQU_CONV.SU = 0.088*CI.FM.HP{1}.GEQU_CONV.Ugs;
CI.FM.HP{1}.GEQU_CONV.UC = CI.FM.HP{1}.GEQU_CONV.Ugs/0.8;
CI.FM.HP{1}.GEQU_CONV.CFL = 0.5;
CI.FM.HP{1}.GEQU_CONV.F_Cutoff = 750;
CI.FM.HP{1}.GEQU_CONV.y_vec = linspace(CI.FM.HP{1}.GEQU_CONV.ra,CI.FM.HP{1}.GEQU_CONV.rb,CI.FM.HP{1}.GEQU_CONV.nb_points);
CI.FM.HP{1}.GEQU_CONV.Lf =  sqrt((CI.FM.HP{1}.GEQU_CONV.Ugs/CI.FM.HP{1}.GEQU_CONV.SU)^2 - 1)*(CI.FM.HP{1}.GEQU_CONV.rb-CI.FM.HP{1}.GEQU_CONV.ra);
CI.FM.HP{1}.GEQU_CONV.dr = 2*CI.FM.HP{1}.GEQU_CONV.Lf/CI.FM.HP{1}.GEQU_CONV.nf_points;
CI.FM.HP{1}.GEQU_CONV.xi_steady = 0:dr:2*Lf;
%CI.FM.HP{1}.GEQU_CONV.uflow =
%CI.FM.HP{1}.GEQU_CONV.xi =

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CI = Fcn_BC(CI)
CI.R1 = 1;
CI.R2 = -1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DET = MyFunDet_NoFlow(s, CI)
R1      = CI.R1;
R2      = CI.R2;
tauP    = CI.TP.tau;
tauM    = CI.TP.tau;
alpha   = CI.TP.alpha;
N       = length(CI.CD.xSp) - 1;
MT      = zeros(2*N);
idf     = CI.CD.indexFlame;
%
for ss = 1 : N - 1
    D = diag([exp(-s*tauP(ss)), exp(s*tauM(ss))]);
    A1 = CI.TPM.B1{ss}*CI.TPM.C1*D;
    A2 = CI.TPM.B2{ss}*CI.TPM.C2;
    MT((ss - 1)*2 + (1 : 2), (ss - 1)*2 + (1 : 2))   = A1;
    MT((ss - 1)*2 + (1 : 2), ss*2 + (1 : 2))         = -A2;
end
%
ep1 = exp(-s*tauP(idf - 1));
ep2 = exp( s*tauM(idf - 1));
ss = idf - 1;
FTF      = CI.FM.n .* exp(-s*CI.FM.td);
MT((ss - 1)*2 + 2, (ss - 1)*2 + 1)   = MT((ss - 1)*2 + 2, (ss - 1)*2 + 1) + alpha .* FTF .* ep1;
MT((ss - 1)*2 + 2, (ss - 1)*2 + 2)   = MT((ss - 1)*2 + 2, (ss - 1)*2 + 2) - alpha .* FTF .* ep2;
% left boundary condition
MT(end - 1, 1 : 2)    = [1, -R1];
MT(end, end - 1 : end)        = [-R2.*exp(-s*tauP(end)), exp( s*tauM(end))];
%
MTsparse = sparse(MT);
DET = det(MTsparse);%./CI.pValWOF;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function CI = Fcn_PreProcessing(CI)
%
CI = Fcn_calculation_Matrix_elements_WOMeanFlow(CI);
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function CI = Fcn_calculation_Matrix_elements_WOMeanFlow(CI)
%
CI.TPM.C1 =  [      1   1;...
    1  -1];
CI.TPM.C2 =  [      1   1;...
    1  -1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(CI.CD.xSp) - 1;
for ss = 1 : N - 1
    [B1,B2] = Fcn_Matrix_calculation_WOMeanFlow(ss, CI);
    CI.TPM.B1{ss}         = B1;               % first step
    CI.TPM.B2{ss}         = B2;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function [B1,B2] = Fcn_Matrix_calculation_WOMeanFlow(ss, CI)
%
Theta   = CI.CD.Theta(ss);
rho1    = CI.TP.rhoMean(ss);
rho2    = CI.TP.rhoMean(ss + 1);
c1      = CI.TP.cMean(ss);
c2      = CI.TP.cMean(ss + 1);
Xi      = Theta*(rho1*c1)/(rho2*c2);
% ----------------------------------
B1(1,1) =   1;
B1(1,2) =   0;
B1(2,1) =   0;
B1(2,2) =   1;
%
B2(1,1) =   1;
B2(1,2) =   0;
B2(2,1) =   0;
B2(2,2) =   Xi;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ev = FcnCalEigVal1(myFun, FreqSp, GrSp, CI, op)
% This function is used to calculate the Eigenvalues of an equation myFun =
% 0;
% FreqSp and GrSp are the frequency samples and growth rate samples
% CI contains the parameters for the function myFun
% op contains the parameters for this solver
% op.prec % precision of the solver, default : 1e-1
% ev is the eigenvalue
if nargin < 5
    op.prec = 1e-2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = optimset('Display','off');        % the calculation results by fsolve are not shown in the workspace

% options = optimoptions('fsolve','FiniteDifferenceType','central', 'Display','off')

for ss  = 1 : length(GrSp)
    for kk  = 1 : length(FreqSp)
        s0  = GrSp(ss) + 1i*2*pi*FreqSp(kk);    % guess value
        F   = @(s)myFun(s, CI);
        [a(ss, kk), ~, flag]   = fsolve(F, s0, options);   % solve equation
        ap(ss, kk)  = floor(a(ss,kk)./op.prec).*op.prec;    % this is used to set the precision
        if flag ~= 1
            ap(ss, kk) = 0;
        end
    end
end
ev = reshape(ap, 1, size(ap, 1)*size(ap, 2));
ev = unique(ev);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the previous processing is still not enough to get rid of the same
% value,such as 999.9 and 1000.1 corresponding to the function `floor'
% such as 995.1 and 994.9 corresponding to the function `round'
cal     = 0;
ev      = sort(ev);
ev_diff = diff(ev);
index_NULL = [];
for kk  = 1 : length(ev_diff)
    if(abs(ev_diff(kk)) < 0.1)
        cal = cal + 1;
        index_NULL(cal) = kk;
    end
end
ev(index_NULL) = []; % this is the EigVal we want
GrLM    = [GrSp(1)    GrSp(end)];
OMLM    = [FreqSp(1)  FreqSp(end)].*2.*pi;
%-------------------------------------
s_null  = [];
cal     = 0;
for ss  = 1 : length(ev)
    if(     real(ev(ss)) < GrLM(1)...
            ||real(ev(ss)) > GrLM(2)...
            ||imag(ev(ss)) < OMLM(1)...
            ||imag(ev(ss)) > OMLM(2)...
            ||imag(ev(ss)) <= 0)
        cal = cal + 1;
        s_null(cal) = ss;
    end
end
ev(s_null) = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
