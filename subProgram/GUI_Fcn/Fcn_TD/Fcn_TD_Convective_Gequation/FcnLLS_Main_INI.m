function [dataGhost, gridGhost, SET, vyPrime_no_stencil] = FcnLLS_Main_INI(vMean, vConv, flameSpeed, dx, dy, ra, rb)
%%
% author: Jingxuan Li
% last updated: 2017-02-26
% initialize the constant data
    SET                     = FcnLLS_INI_CstData(vMean, vConv, flameSpeed, dx, dy, ra, rb);
    % grid and signed function initialization
    [dataGhost, gridGhost]  = FcnLLS_INI_Data_Grid(SET.step, SET.RA, SET.RB, SET.sL0, SET.vMean, SET.vConv);
    pause(1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % steady information
    %     % search the indices of the local data
%     [LLSid, dataLLS]        = FcnLLS_INI(dataGhost, gridGhost);
%     FF0                     = FcnLLS_FrontTrack_Summary_Steady_Flame(tNow, dataLLS, dataGhost, LLSid, gridGhost, SET);
    % initialise vy prime
    vyPrime_no_stencil      = Fcn_INI_VFresh(gridGhost);
%
end
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  breaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SET = FcnLLS_INI_CstData(vMean, vConv, flameSpeed, dx, dy, ra, rb)
%%
% author: Jingxuan Li
% last updated: 2017-02-26
% This function is used to set the constant data for the calculation
% These data cannot be changed, once the calculation has been started
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % basic flame speed
    SET.sL0     = flameSpeed;
    % basic flow velocity
    SET.vMean   = vMean;
    %AO: convective model speed
    SET.vConv   = vConv;
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % spatial information
    % spatial step
    SET.step    = [dx dy];
    % A%
    % radius of the internal burner
    SET.RA       = ra;
    % radius of external burner
    SET.RB       = rb;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % time information
    % time step
%    SET.dt      = dt;     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % reinitialisation information
    SET.ReINI.NStep = 4;
 %   SET.ReINI.dt    = 0.5*SET.dt;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % saving name
    SET.SaveName = ['case_1'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function vyPrime_no_stencil = Fcn_INI_VFresh(gridGhost)
%%
% author: Jingxuan Li
% last updated: 2017-02-26
    ns                  = gridGhost.stencil;
    %
    y_no_stencil        = gridGhost.vs{2}(ns + 1 : end - ns);
    vyPrime_no_stencil  = 0.*y_no_stencil;
    sz = size(vyPrime_no_stencil);
    if sz(2) < sz(1)
        vyPrime_no_stencil = vyPrime_no_stencil.'; 
    else
        vyPrime_no_stencil = vyPrime_no_stencil;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  breaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataGhost, gGhost] = FcnLLS_INI_Data_Grid(step, RA, RB, sLMean, vMean, vConv)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
if nargin < 1
    step = 5*[1 1].*1e-5;
end
if nargin < 2
    RA = 0; 
end
if nargin < 3
   RB = 10e-3; 
end
if nargin < 4
   sLMean = 0.368; 
end
if nargin < 5
   vMean = 2.12; 
end
if nargin < 6
    vConv = vMean;
end
%
    [data, dataGhost, g, gGhost] = FcnLLS_INI_Data_Grid_INI(step, RA, RB, sLMean, vMean);
    %
%     w  = abs(gGhost.xs{1}(end) - gGhost.xs{1}(1));
%     h  = abs(gGhost.xs{2}(end) - gGhost.xs{2}(1));
%     asp = h./w;
%     h = figure;
%     scrsz = get(0,'ScreenSize');
%     set(h,'Position',[scrsz(4).*(1/8) scrsz(4).*(1/20) scrsz(3)*2/5 scrsz(4).*(4/5)])
%     hAxes(1)=axes('Unit','pixels','position',[100 100 400 400.*asp./3]);
%     contour(gGhost.xs{1},gGhost.xs{2},dataGhost, 50, 'linewidth',1);
%     xlabel('x [m]')
%     ylabel('y [m]')
%     grid on
%     title('Initial G')
%     colorbar
    assignin('base','g',g);
    assignin('base','gGhost',gGhost);
    assignin('base','data',data);
    assignin('base','dataGhost',dataGhost);
    %
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
function [data, dataGhost, g, gGhost] = FcnLLS_INI_Data_Grid_INI(step, RA, RB, sL, vMean)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
    stencil     = 4;
    xLim        = (RB+2*stencil*step(1))*[-1 0];
    Lf          = sqrt(vMean^2/sL^2-1)*(RB-RA);
    yLim        = [0 Lf]*1.75;
    xLimGhost   = xLim + stencil*step(1).*[-1 1];
    yLimGhost   = yLim + stencil*step(2).*[-1 1];
    g           = GridMaking(xLim, yLim, step);
    gGhost      = GridMaking(xLimGhost, yLimGhost, step);
    %%%%%%%%%%%%%%%%%%%
    %% AO: Use this for a conical flame
%     data        = FcnSignedDistINI(g, RB, sL, vMean);
%     dataGhost   = FcnSignedDistINI(gGhost, RB, sL, vMean);
%     g           = FcnLocateFlameBase(g, RB);
%     gGhost      = FcnLocateFlameBase(gGhost, RB);
    %% AO: Use this for a V-flame
    data        = FcnSignedDistINI_VFlame(g, RA, RB, sL, vMean);
    dataGhost   = FcnSignedDistINI_VFlame(gGhost, RA, RB, sL, vMean);
    g           = FcnLocateFlameBase(g, RA);
    gGhost      = FcnLocateFlameBase(gGhost, RA);
    gGhost.stencil  = stencil;
    g.stencil   = stencil;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function data = FcnSignedDistINI(g, R, sL, vMean)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
    alpha   = asin(sL/vMean);
    H       = R/tan(alpha);
    h       = R*tan(alpha);
    k1      = 1/tan(alpha);
    k2      = -1/k1;
    xgrid = g.xs{1};
    ygrid = g.xs{2};
    data = zeros(size(g.xs{1}));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ValLine3L = g.xs{2} - k1.*g.xs{1} - H;
    ValLine1L = g.xs{2} - k2.*g.xs{1} - H;
    ValLine2L = g.xs{2} - k2.*g.xs{1} + h;
    ValLine3R = g.xs{2} + k1.*g.xs{1} - H;
    ValLine1R = g.xs{2} + k2.*g.xs{1} - H;
    ValLine2R = g.xs{2} + k2.*g.xs{1} + h;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dist1 = ((g.xs{1} - 0).^2 + (g.xs{2} - H).^2).^0.5;         % distance to the flame tip     
    Dist2 = ((g.xs{1} + R).^2 + (g.xs{2} - 0).^2).^0.5;         % distance to the flame left base
    Dist3 = ((g.xs{1} - R).^2 + (g.xs{2} - 0).^2).^0.5;         % distance to the flame right base
    normDen = (1 + k1.^2).^0.5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ss = 1 : size(data, 1)
        for kk = 1 : size(data, 2)
            if     xgrid(ss, kk) <  -R && ValLine2L(ss, kk) <= 0
                data(ss, kk) =  Dist2(ss,kk);
            elseif xgrid(ss, kk) >= -R && xgrid(ss, kk) <= 0 && ValLine2L(ss, kk) <= 0
                data(ss, kk) = -Dist2(ss,kk);
            elseif xgrid(ss, kk) <= 0 && ValLine2L(ss, kk) >  0 && ValLine1L(ss, kk) <= 0
                data(ss, kk) = ValLine3L(ss,kk)./normDen;
            elseif xgrid(ss, kk) <= 0 && ValLine1L(ss, kk) >  0
                data(ss, kk) = Dist1(ss,kk);
            elseif xgrid(ss, kk) <= R && xgrid(ss, kk) >  0  && ValLine2R(ss, kk) <= 0
                data(ss, kk) = -Dist3(ss,kk);
            elseif xgrid(ss, kk) >  R  && ValLine2R(ss, kk) <= 0
                data(ss, kk) =  Dist3(ss,kk);
            elseif xgrid(ss, kk) >  0 && ValLine2R(ss, kk) >  0 && ValLine1R(ss, kk) <= 0
                data(ss, kk) = ValLine3R(ss,kk)./normDen;
            elseif xgrid(ss, kk) >  0 && ValLine1R(ss, kk) >  0
                data(ss, kk) = Dist1(ss,kk);
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%


function data = FcnSignedDistINI_VFlame(g, RA, R, sL, vMean)
%%
% author: Alessandro Orchini
% last updated: 2017-02-22
    H       = sqrt(vMean^2/sL^2-1)*(R-RA);
%    alpha   = atan((R-RA)/H);
    k1      = H/(-R+RA);
    k2      = -1/k1;
    xgrid   = g.xs{1};
    data    = zeros(size(g.xs{1}));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ValLine3L = g.xs{2} - k1.*g.xs{1} - H + k1*(-R);
    ValLine1L = g.xs{2} - k2.*g.xs{1} - H + k2*(-R);
    ValLine2L = g.xs{2} - k2.*g.xs{1}     + k2*(-RA);
    ValLine3R = g.xs{2} + k1.*g.xs{1} - H + (-k1)*R;
    ValLine1R = g.xs{2} + k2.*g.xs{1} - H + (-k2)*R;
    ValLine2R = g.xs{2} + k2.*g.xs{1}     + (-k2)*RA;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dist1a = ((g.xs{1} + R).^2 + (g.xs{2} - H).^2).^0.5;         % distance to the flame tip left     
    Dist1b = ((g.xs{1} - R).^2 + (g.xs{2} - H).^2).^0.5;         % distance to the flame tip right     
    Dist2 = ((g.xs{1} + RA).^2 + (g.xs{2} - 0).^2).^0.5;         % distance to the flame left base
    Dist3 = ((g.xs{1} - RA).^2 + (g.xs{2} - 0).^2).^0.5;         % distance to the flame right base
    normDen = (1 + k1.^2).^0.5;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ss = 1 : size(data, 1)
        for kk = 1 : size(data, 2)
                %% Left Part
            if     xgrid(ss,kk) <=0 && xgrid(ss, kk) >-RA && ValLine2L(ss, kk) <= 0
                %1
                data(ss, kk) =  Dist2(ss,kk);
            elseif xgrid(ss, kk) <= -RA && ValLine2L(ss, kk) <= 0
                %2
                data(ss, kk) = -Dist2(ss,kk);
            elseif xgrid(ss, kk) <= 0 && ValLine2L(ss, kk) >  0 
                %&& ValLine1L(ss, kk) <= 0
                %3
                data(ss, kk) = ValLine3L(ss,kk)./normDen;
%             elseif xgrid(ss, kk) <= 0 && ValLine1L(ss, kk) >  0 && ValLine3L(ss,kk) > 0 
%                 %4
%                 data(ss, kk) = Dist1a(ss,kk);
%             elseif xgrid(ss, kk) <= 0 && ValLine1L(ss, kk) >  0 && ValLine3L(ss,kk) <= 0 
%                 %4bis
%                 data(ss, kk) = -Dist1a(ss,kk);
                %% Right part
            elseif xgrid(ss, kk) >  RA  && ValLine2R(ss, kk) <= 0
                %5
                data(ss, kk) = -Dist3(ss,kk);
            elseif xgrid(ss, kk)>0 && xgrid(ss, kk) < RA  && ValLine2R(ss, kk) <= 0
                %6
                data(ss, kk) =  Dist3(ss,kk);
            elseif xgrid(ss, kk) >  0 && ValLine2R(ss, kk) >  0 
                %&& ValLine1R(ss, kk) <= 0
                %7
                data(ss, kk) = ValLine3R(ss,kk)./normDen;
%             elseif xgrid(ss, kk) >  0 && ValLine1R(ss, kk) >  0 && ValLine3L(ss,kk) > 0
%                 %8
%                 data(ss, kk) = Dist1b(ss,kk);
%             elseif xgrid(ss, kk) >  0 && ValLine1R(ss, kk) >  0 && ValLine3L(ss,kk) <= 0
%                 %8bis
%                 data(ss, kk) = -Dist1b(ss,kk);
            end
        end
    end
end

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function g = FcnLocateFlameBase(g, R)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
    xgrid = g.xs{1};
    ygrid = g.xs{2};
    xSp = xgrid(1,:);
    [~,g.indexBaseX(1)] = min(abs(xSp - (-R)));
    ySp = ygrid(:,1);
    [~,g.indexBaseY(1)] = min(abs(ySp));
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function gridOut = GridMaking(xLim, yLim, step)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
%
    gridOut.dim     = 2;
    gridOut.min     = [xLim(1); yLim(1)];
    gridOut.max     = [xLim(2); yLim(2)];
    gridOut.dx      = [step(1); step(2)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gridOut.vs = cell(gridOut.dim, 1);
    for i = 1 : gridOut.dim
      gridOut.vs{i} = (gridOut.min(i) : gridOut.dx(i) : gridOut.max(i))';
    end
    %%%%%%%%%%%%%%%%%%
    gridOut.N = zeros(gridOut.dim, 1);
    for i = 1 : gridOut.dim
      gridOut.N(gridOut.dim + 1 - i) = length(gridOut.vs{i});
    end
    %%%%%%%%%%%%%%%%%%
    gridOut.xs = cell(gridOut.dim, 1);
    if(gridOut.dim > 1)
      [gridOut.xs{:}] = meshgrid(gridOut.vs{:});
    end
    %%%%%%%%%%%%%%%%%%
    gridOut.axis = zeros(1, 2 * gridOut.dim);
    for i = 1 : gridOut.dim
        gridOut.axis(2 * i - 1 : 2 * i) = [ gridOut.min(i), gridOut.max(i) ];
    end
    if(gridOut.dim == 1)
        gridOut.shape = [ gridOut.N, 1 ];
    else
        gridOut.shape = gridOut.N';
    end  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
