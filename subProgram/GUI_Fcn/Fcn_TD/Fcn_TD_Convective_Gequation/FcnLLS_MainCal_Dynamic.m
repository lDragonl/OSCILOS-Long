function [dataGhost, FF, vyPrime_no_stencil, tNow] = ...
    FcnLLS_MainCal_Dynamic(tNow, SET, dataGhost, gridGhost,...
    vyPrime_no_stencil, vyPrime_inlet, isReINI_Global)
%%
% author: Jingxuan Li
% last updated: 2017-02-26
if nargin < 7
    isReINI_Global = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% search the indices of the local data
[LLSid, dataLLS]        = FcnLLS_INI(dataGhost, gridGhost);
% main calculation step
tSpan                   = [tNow, tNow + SET.dt];
%
[vyPrime_no_stencil]    = Fcn_Current_vy_Prime_No_Stencil(vyPrime_no_stencil, vyPrime_inlet, tNow, gridGhost.dx(2), SET.dt, SET.vConv);
%
[tNow, dataLLS, Info]   = FcnTVDRunge_Kutta3(tSpan, dataLLS, dataGhost, gridGhost, SET, LLSid, vyPrime_no_stencil);
% reinitialisation
dataLLS                 = FcnLLS_ReINI(dataLLS, dataGhost, gridGhost, LLSid, SET.ReINI.NStep, SET.ReINI.dt);
% update the global data
dataGhost               = FcnLLS2dataGhost(dataLLS, dataGhost, LLSid, gridGhost);
%
% to avoide the divergence, a global reinitialisation is needed
if isReINI_Global == 1
    dataGhost           = FcnLLS_ReINI_Global(dataGhost, gridGhost, SET.ReINI.NStep, SET.ReINI.dt);
end
%
FF                      = FcnLLS_FrontTrack_Summary(tNow, dataLLS, dataGhost, LLSid, gridGhost, Info, SET);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%c%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  breaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vyPrime_no_stencil, tNow] = Fcn_Current_vy_Prime_No_Stencil(vyPrime_no_stencil, vyPrime_inlet, tNow, dy, dt, vMean)
%%
% author: Jingxuan Li
% last updated: 2017-02-26
[vyPrime_no_stencil, tNow]  = Fcn_1D_Convection_Combine(vyPrime_no_stencil, tNow, dy, dt, vMean);
vyPrime_no_stencil(1)       = vyPrime_inlet; % update the inlet element
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vNew, tNew] = Fcn_1D_Convection_Combine(vOld, tOld, dx, dt, vMean)
%%
% author: Jingxuan Li
% last updated: 2017-02-26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Solving 1-D wave equation with WENO5 (MOL-WENO5-LF)
%                 dv/dt + vMean*dv/dx = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fydot           = @(y, dx, vMean)Fcn_WENO5_Advection(y, dx, vMean);
tSpan           = tOld + [0, dt];
[tNew, vNew]    = Fcn_TVD_RK3(tSpan, vOld, Fydot, dx, vMean);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, y] = Fcn_TVD_RK3(tSpan, y0, Fydot, dx, vMean)
%%
% author: Jingxuan Li
% last updated: 2017-02-24
% vMean is the mean axial velocity
% Fydot is the function to calculate ydot
%
t       = tSpan(1);
y       = y0;
dt      = diff(tSpan);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydot    = feval(Fydot, y, dx, vMean);
%
t1      = t + dt;
y1      = y + dt*ydot;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydot    = feval(Fydot, y1, dx, vMean);
% Take the second substep.
t2      = t1 + dt;
y2      = y1 + dt*ydot;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tHalf   = 0.25*(3*t + t2);
yHalf   = 0.25*(3*y + y2);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydot    = feval(Fydot, yHalf, dx, vMean);
tThreeHalf  = tHalf + dt;
yThreeHalf  = yHalf + dt*ydot;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t       = (1/3)*(t + 2*tThreeHalf);
y       = (1/3)*(y + 2*yThreeHalf);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = Fcn_WENO5_Advection(y0, dx, vMean)
%%
% author: Jingxuan Li
% last updated: 2017-02-24
%
sz = size(y0);
if sz(2) < sz(1)
    y = y0.';
else
    y = y0;
end
N = length(y);
ym          = circshift(y, [0 1]);
dydx0       = (y - ym)./dx;
D1(1, :)    = circshift(dydx0, [0 2]);
D1(2, :)    = circshift(dydx0, [0 1]);
D1(3, :)    = dydx0;
D1(4, :)    = circshift(dydx0, [0 -1]);
D1(5, :)    = circshift(dydx0, [0 -2]);
%
dydx        = FcnUpWindFirstWENO5_one_direction(D1);
%%%%%%%%%%%%%%%%%%%%%%%
dydx(1)     = (y(2) - y(1)) ./ dx;
dydx(2)     = Fcn_dydx_Center(y(1 : 3), dx);
dydx(3)     = Fcn_dydx_Center(y(2 : 4), dx);
dydx(N - 1) = Fcn_dydx_Up(y(N - 3 : N - 1), dx);
dydx(N)     = Fcn_dydx_Up(y(N - 2 : N), dx);
%%%%%%%%%%%%%%%%%%%%%%%
dydt        = -vMean .* dydx;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydx = Fcn_dydx_Center(y, dx)
%%
% author: Jingxuan Li
% last updated: 2017-02-24
% central second order differential
dydx = (y(3) - y(1)) ./ dx ./ 2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydx = Fcn_dydx_Up(y, dx)
%%
% author: Jingxuan Li
% last updated: 2017-02-24
% upwind second order differential
dydx = (1.5*y(3) - 2*y(2) + 0.5*y(1)) ./ dx;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  breaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FFNew = FcnLLS_FrontTrack_Summary(tNow, dataLLS, dataGhost, LLSid, gridGhost, Info,  SET)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% This function is used to save the flame front information
[dataFF, colFF, rowFF, idxLLSFF] = FcnLLS_FrontTrack(dataLLS, dataGhost, LLSid, gridGhost);
%
FF.tNow         = tNow;
% information
FF.curvature    = Info.curvature(idxLLSFF);
FF.velocity{1}  = Info.velocity{1}(idxLLSFF);
FF.velocity{2}  = Info.velocity{2}(idxLLSFF);
FF.sL           = Info.sL(idxLLSFF);
FF.cof          = Info.cof(idxLLSFF);
%
FF.dataFF       = dataFF;
FF.colFF        = colFF;
FF.rowFF        = rowFF;
FF.idxLLSFF     = idxLLSFF;
%idxFF           = sub2ind(size(dataGhost), rowFF, colFF);
idxFF           = size(dataGhost,1)*(colFF*1)+rowFF;
FF.xFF          = gridGhost.xs{1}(idxFF);
FF.yFF          = gridGhost.xs{2}(idxFF);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
FFNew           = FcnFFFrontLocate(FF, SET, dataGhost, gridGhost);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataFF, colFF, rowFF, idxLLSFF] = FcnLLS_FrontTrack(dataLLS, dataGhost, LLSid, gridGhost)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% This function is used to track the flame front
%
alpha           = 2*max(gridGhost.dx);
% find the indices of the data being suspected to be the elements of Gamma
% idxLLSFF        = find( abs(dataLLS(1 : LLSid.NMargin(2))) <=  alpha );
% AO: change 1
idxLLSFF        = find( abs(dataLLS(1 : LLSid.NMargin(3))) <=  alpha );
nFF             = length(idxLLSFF);
rowFF           = LLSid.rowLLS(idxLLSFF);
colFF           = LLSid.colLLS(idxLLSFF);
%idxFF           = sub2ind(size(dataGhost), rowFF, colFF);
sz              = size(dataGhost,1);
idxFF           = sz*(colFF-1)+rowFF;
dataFF          = dataGhost(idxFF)';   % The data being suspected to be the elements of Gamma
dataFFN         = zeros(4, nFF);
%
isFFN           = zeros(4, nFF);
%
% idxFFN(1,:)     = sub2ind(size(dataGhost), rowFF - 1, colFF);
% idxFFN(2,:)     = sub2ind(size(dataGhost), rowFF + 1, colFF);
% idxFFN(3,:)     = sub2ind(size(dataGhost), rowFF, colFF - 1);
% idxFFN(4,:)     = sub2ind(size(dataGhost), rowFF, colFF + 1);
idxFFN(1,:)     = rowFF - 1 + (colFF     -1)*sz;
idxFFN(2,:)     = rowFF + 1 + (colFF     -1)*sz;
idxFFN(3,:)     = rowFF     + (colFF - 1 -1)*sz;
idxFFN(4,:)     = rowFF     + (colFF + 1 -1)*sz;
%
for kk = 1 : 4
    dataFFN(kk,:)    = dataGhost( idxFFN(kk,:) );
    isFFN(kk,:)      = ( dataFFN(kk,:).*dataFF <= 0 );
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
isFFNsum            = sum(isFFN, 1);
% find the cell not belonging to Gamma
isNoFFN             = find(isFFNsum == 0);
%
idxLLSFF(isNoFFN)   = [];
% nFF                 = nFF - length(isNoFFN);
rowFF(isNoFFN)      = [];
colFF(isNoFFN)      = [];
idxFF(isNoFFN)      = [];
dataFF(isNoFFN)     = [];
%
end
%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------

function FFNew = FcnFFFrontLocate(FF, SET, dataGhost, gridGhost)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
FFNew.tNow = FF.tNow;
%dx      = min(SET.step);
%xFF     = FF.xFF;
%yFF     = FF.yFF;
%dataFF  = FF.dataFF;
%[~, FFNew.indicesYNew]  = sortrows([yFF,-xFF]);
%FFNew.xFF               = xFF(FFNew.indicesYNew);
%FFNew.yFF               = yFF(FFNew.indicesYNew);
%FFNew.dataFF            = dataFF(FFNew.indicesYNew);
%     [FFNew.xFFInterp, FFNew.yFFInterp]...
%         = FcnFFFrontInterp( FFNew.xFF, FFNew.yFF, FFNew.dataFF, dx);
%     FFNew.A = FcnFFSurfaceAreaNew(FFNew.xFFInterp, FFNew.yFFInterp);
[FFNew.xFFInterp, FFNew.yFFInterp] = FcnFFFrontInterpNew( dataGhost, gridGhost);
FFNew.A = 0 ;
for ii = 1:length(FFNew.xFFInterp)
%     [~, FFNew.indicesYNew]  = sortrows([FFNew.yFFInterp{ii}',-FFNew.xFFInterp{ii}']);
%     FFNew.xFFInterp{ii}     = FFNew.xFFInterp{ii}(FFNew.indicesYNew);
%     FFNew.yFFInterp{ii}     = FFNew.yFFInterp{ii}(FFNew.indicesYNew);
    FFNew.A                 = FFNew.A + FcnFFSurfaceAreaNew(FFNew.xFFInterp{ii}, FFNew.yFFInterp{ii});
end

%[surfArea, deltau] = delta_Int( dataGhost, gridGhost, FF );
%FFNew.Anew = surfArea;
%figure(77)
%contourf(gridGhost.xs{1},gridGhost.xs{2},deltau)
%drawnow
end

function [xFFInterp, yFFInterp] = FcnFFFrontInterpNew( dataGhost, gridGhost)

a = contourcs(gridGhost.xs{1}(1,1+3*gridGhost.stencil:end-gridGhost.stencil),gridGhost.xs{2}(1+gridGhost.stencil:end-gridGhost.stencil,1),dataGhost(1+gridGhost.stencil:end-gridGhost.stencil,1+3*gridGhost.stencil:end-gridGhost.stencil),[0 0]);

for ii = 1:length(a)
    a(ii).X(a(ii).Y==0) = -gridGhost.RA;
    xFFInterp{ii} = a(ii).X;
    yFFInterp{ii} = a(ii).Y;
end




end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xFFFinal, yFFFinal] = FcnFFFrontInterp(xFF, yFF, dataFF, dx)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% currently, we donot account for the change of flame speed sL, because it
% is not easy to interpolate it.
N   = length(yFF);          % number of candidates
%%%%%%%%%%%%%%%%%%%%%% first step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first step: interpolate the points on zero level sets from y
% direction. The objective of this step is to find the points of
% intersection for the zero level sets and the x grids
ss = 1;    % The index of point with the smallest value of x for a given y
kk = 1;    % The index of envisaged y
while ss <= N
    [~, Ind] = find(yFF == yFF(ss));
    ny = length(Ind);  % number of points with the same value of x for a given y
    if ny > 1
        % zero level sets cut through this x grid
        yFFInterp1(kk)   = yFF(ss);
        % we currently use the first order approximation. Linear
        % interpolation is used.
        xFFInterp1(kk)   = interp1(  dataFF(ss : ss + ny - 1),...
            xFF(ss : ss + ny - 1), 0, 'linear','extrap');
        indexY(kk)  = ss;    % index of the horizontal grid
        numY(kk)    = ny;    % number of points with the same value of x
        ss = ss + ny;
        kk = kk + 1;
    else
        %       % zero level sets pass though this point
        yFFInterp1(kk)      = yFF(ss);
        xFFInterp1(kk)      = xFF(ss);
        indexY(kk)          = ss;    % index of the horizontal grid
        numY(kk)            = ny;    % number of points with the same value of x
        ss = ss + ny;
        kk = kk + 1;
        
    end
end
% at the end of the flame, set the x position to be 0
% AO: this is not needed for the V-flame
%xFFInterp1(end) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%% Second step %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The second step is to find the points of intersection for the zero level
% sets and the y grids
xFFFinal    = xFFInterp1;
yFFFinal    = yFFInterp1;
Nx          = xFFInterp1./dx;
NxLeft      = floor(Nx);            % Index of the left grid of the cut point
nxPool      = abs(diff(NxLeft));    % number of points at the vertical grid between each neighbour horizontal grids
id = 0;                             % index of final element
N1 = length(xFFInterp1);
for ss = 1 : N1 - 1
    id = id + 1;
    xFFFinal(id)    = xFFInterp1(ss);
    yFFFinal(id)    = yFFInterp1(ss);
    if nxPool(ss)   ~= 0
        for kk = 1 : nxPool(ss)
            if NxLeft(ss + 1) > NxLeft(ss)
                xFFFinal(id + kk) = (NxLeft(ss) + kk).*dx;
            else
                xFFFinal(id + kk) = (NxLeft(ss) + 1 - kk).*dx;
            end
            %           disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
            [~, Ind1]   = min(abs(xFF(indexY(ss) : indexY(ss) + numY(ss) - 1) - xFFFinal(id + kk)));
            y1(1)       = yFF(indexY(ss) - 1 + Ind1);
            dataFF1(1)  = dataFF(indexY(ss) - 1 + Ind1);
            [~, Ind2]   = min(abs(xFF(indexY(ss + 1) : indexY(ss + 1) + numY(ss + 1) - 1) - xFFFinal(id + kk)));
            y1(2)       = yFF(indexY(ss + 1) - 1 + Ind2);
            dataFF1(2)  = dataFF(indexY(ss + 1) - 1 + Ind2);
            yFFFinal(id + kk) = interp1(  dataFF1, y1, 0, 'linear','extrap');
        end
    end
    id = id + nxPool(ss);
end
N = length(yFFFinal);
%     try % treatment of the flame end
%         xFlameEnd   = xFF(indexY(end) : end);
%         yFlameEnd   = yFF(indexY(end) : end);
%         dataFlameEnd = dataFF(indexY(end) : end);
%         %
%         y1(1)       = yFF(indexY(end) - 1);
%         dataFF1(1)  = dataFF(indexY(end) - 1);
%         [~, Ind2]   = min(abs(xFlameEnd));
%         y1(2)       = yFlameEnd(Ind2);
%         dataFF1(2)  = dataFlameEnd(Ind2);
%         xFFFinal(N + 1) = 0;
%         yFFFinal(N + 1) = interp1(  dataFF1, y1, 0, 'linear','extrap');
%     end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = FcnFFSurfaceAreaNew(xFF, yFF)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% calculate Dx
% calculate ds
Dx = abs((xFF(1:end-1) + xFF(2:end))./2);
dx = abs(diff(xFF));
dy = abs(diff(yFF));
indX = (dx > dy);
indY = (dx <= dy);
ds = indY.*dy.*(1 + (dx./(dy + eps)).^2).^0.5 + indX.*dx.*(1 + (dy./(dx + eps)).^2).^0.5;
A = 2*pi*sum(Dx.*ds);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  breaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [curvatureMean, gradMagMean] = FcnLLS_curvatureMean(dataGhost, LLSid, gridGhost)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
sizeSquare  = 1;
[Irow,Jcol] = ind2sub(sizeSquare*[1, 1], (1 : sizeSquare^2));
idxCol      = Jcol(:) - (sizeSquare + 1)/2;
idxRow      = Irow(:) - (sizeSquare + 1)/2;
%
phix        = zeros(sizeSquare^2, length(LLSid.rowLLS));
phiy        = phix;
phixy       = phix;
phixx       = phix;
phiyy       = phix;
for kk = 1 : length(idxCol)
    [derivFirst, derivSecond] = FcnLLS_curvature_5Items(dataGhost, LLSid, gridGhost, idxCol(kk), idxRow(kk));
    phix(kk,:)  = derivFirst{2};
    phiy(kk,:)  = derivFirst{1};
    phixy(kk,:) = derivSecond{1,2};
    phixx(kk,:) = derivSecond{2,2};
    phiyy(kk,:) = derivSecond{1,1};
end
phixMean        = mean(phix, 1);
phiyMean        = mean(phiy, 1);
phixyMean       = mean(phixy, 1);
phixxMean       = mean(phixx, 1);
phiyyMean       = mean(phiyy, 1);
phixPowerMean   = phixMean.^2;
phiyPowerMean   = phiyMean.^2;
phiMixMean      = phixMean.*phiyMean;
gradMagMean     = (phixPowerMean + phiyPowerMean).^0.5;
%---------------------------------------------------------------------------
curvatureMean   = zeros(size(gradMagMean));
curvatureMean   = curvatureMean + phixxMean .* phiyPowerMean...
    - 2*phixyMean .* phiMixMean ...
    + phiyyMean .* phixPowerMean;
nonzero         = find(gradMagMean > 0);
curvatureMean(nonzero)  = curvatureMean(nonzero) ./ (gradMagMean(nonzero)).^3;
stencil = gridGhost.stencil;
% The curvature of elements near the base of the flame and near the periodic
% boundary are considered to be zero
curvatureMean( LLSid.rowLLS < stencil + 3 ) = 0;
curvatureMean( LLSid.colLLS > gridGhost.N(2) - stencil - 3 ) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [derivFirst, derivSecond] = FcnLLS_curvature_5Items(dataGhost, LLSid, gridGhost, idxCol, idxRow)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% construct the matrix for the curvature calculation
MtrCurv                     = FcnLLS_curvature_Mtr(dataGhost, LLSid, idxCol, idxRow);
% calculate the first and second partial derivative
[derivFirst, derivSecond]   = FcnLLS_curvature_HessianSecond(MtrCurv, gridGhost);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [derivFirst, derivSecond] = FcnLLS_curvature_HessianSecond(Mtr, gridGhost)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% Calculate the first order centered difference approximation
% Calculate the second order centered difference approximation to the Hessian
% matrix (the second order mixed spatial derivative of the data).
%
dim         = 2;
derivFirst  = cell(dim,1);
derivSecond = cell(dim,dim);
dxInv       = 1./gridGhost.dx;
dxInvPower  = dxInv.^2;
dxInvMix    = dxInv(1).*dxInv(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the first order spatial derivative of the data
derivFirst{2}       = 0.5 .* dxInv(2) .* ( Mtr(8, :) - Mtr(2, :) );
derivFirst{1}       = 0.5 .* dxInv(1) .* ( Mtr(6, :) - Mtr(4, :) );
% calculate the second order mixed spatial derivative of the data
derivSecond{2, 2}   = dxInvPower(2) .* ( Mtr(8, :) - 2*Mtr(5, :) + Mtr(2, :) );
derivSecond{1, 1}   = dxInvPower(1) .* ( Mtr(6, :) - 2*Mtr(5, :) + Mtr(4, :) );
derivSecond{1, 2}   = 0.25 .* dxInvMix .* ( Mtr(9, :) + Mtr(1, :) - Mtr(3, :) - Mtr(7, :) );
derivSecond{2, 1}   = derivSecond{1, 2};
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MtrCurv = FcnLLS_curvature_Mtr(dataGhost, LLSid, idxCol, idxRow)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% This function is used to contruct the matrix used in the curvature
% calculation
%
MtrCurv     = zeros(9, length(LLSid.rowLLS));
indexCol    = [-1, -1, -1,  0,  0,  0,  1,  1,  1];
indexRow    = [-1,  0,  1, -1,  0,  1, -1,  0,  1];
%
rowLLSGhost     = LLSid.rowLLS;
colLLSGhost     = LLSid.colLLS;
sz = size(dataGhost,1);
for kk = 1 : 9
%     idx             = sub2ind(size(dataGhost),  rowLLSGhost + indexRow(kk) + idxRow,...
%         colLLSGhost + indexCol(kk) + idxCol);
    idx             = rowLLSGhost + indexRow(kk) + idxRow + ...
                     (colLLSGhost + indexCol(kk) + idxCol - 1)*sz;
    MtrCurv(kk,:)   = dataGhost(idx);
end
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  breaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, y, Info] = FcnTVDRunge_Kutta3(tSpan, y0, dataGhost, gridGhost, SET, LLSid, vyPrime_no_stencil)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
t       = tSpan(1);
y       = y0;
deltaT  = diff(tSpan);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydot    = FcnLLS_ydot(y, t, dataGhost, gridGhost, SET, LLSid, vyPrime_no_stencil);
%
t1      = t + deltaT;
y1      = y + deltaT*ydot;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydot    = FcnLLS_ydot(y1, t1, dataGhost, gridGhost, SET, LLSid, vyPrime_no_stencil);
t2      = t1 + deltaT;
y2      = y1 + deltaT*ydot;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tHalf   = 0.25*(3*t + t2);
yHalf   = 0.25*(3*y + y2);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ydot, Info]  = FcnLLS_ydot(yHalf, tHalf, dataGhost, gridGhost, SET, LLSid, vyPrime_no_stencil);
tThreeHalf  = tHalf + deltaT;
yThreeHalf  = yHalf + deltaT*ydot;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t       = (1/3)*(t + 2*tThreeHalf);
y       = (1/3)*(y + 2*yThreeHalf);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ydot, Info] = FcnLLS_ydot(dataLLS, tNow, dataGhost, gridGhost, SET, LLSid, vyPrime_no_stencil)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% update the dataGhost
dataGhost       = FcnLLS2dataGhost(dataLLS, dataGhost, LLSid, gridGhost);
Mtr             = FcnLLS_MtrWENO(dataGhost, LLSid);
% calculate the spatial derivative using WENO5 approach
[derivL,derivR] = FcnUpWindFirstWENO5(Mtr, gridGhost);
cof             = FcnLLS_cof(dataLLS, LLSid);
%    AO: these are not needed when Mark = 0;
%[curvature, gradMag] = FcnLLS_curvatureMean(dataGhost, LLSid, gridGhost);
curvature = zeros(1,length(LLSid.rowLLS));
gradMag = ones(1,length(LLSid.rowLLS));
[velocity, sL]  = FcnTermVelocity_dynamic(gridGhost, tNow, SET, LLSid, cof, vyPrime_no_stencil);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydot            = FcnTermCV_sL(derivL, derivR, velocity, sL, curvature, gradMag, SET);
Info.curvature  = curvature;
Info.velocity   = velocity;
Info.sL         = sL;
Info.cof        = cof;
end
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cof = FcnLLS_cof(dataLLS, LLSid)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% cut off function: COF
% cof is from the reference Peng et al (1999)
N           = length(LLSid.rowLLS);
% change to 1 x N style
cof         = ones(1, N);
% range2      = [ LLSid.NMargin(2) + 1: LLSid.NMargin(1) ];
% cof(range2) = (abs(dataLLS(range2)) - LLSid.gamma(1)).^2.*...
%     (2*abs(dataLLS(range2)) + LLSid.gamma(1) - 3*LLSid.gamma(2))./...
%     (LLSid.gamma(1) -   LLSid.gamma(2)).^3;
% cof(LLSid.NMargin(1) + 1 : N) = 0;
%AO: change 2
range2      = [ LLSid.NMargin(3) + 1: LLSid.NMargin(2) ];
cof(range2) = (abs(dataLLS(range2)) - LLSid.gamma(2)).^2.*...
    (2*abs(dataLLS(range2)) + LLSid.gamma(2) - 3*LLSid.gamma(3))./...
    (LLSid.gamma(2) -   LLSid.gamma(3)).^3;
cof(LLSid.NMargin(1) + 1 : N) = 0;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Velocity, sL] = FcnTermVelocity_dynamic(gridGhost, tNow, SET, LLSid, cof, vyPrime_no_stencil)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
Velocity    = FcnTermVelocityConvection_dynamic(gridGhost, tNow, SET, LLSid, cof, vyPrime_no_stencil);
% The flame speed may depend on the curvature of the flame front.
% We currently set it to be constant value
sL          = FcnTermsL(LLSid, SET, cof);
%
end
%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function sL = FcnTermsL(LLSid, SET, cof)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
%------------------------------
sL0 = SET.sL0;
% cof is the cut-off function
sL  = cof.*sL0;
%
end
%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
%
function [Velocity] = FcnTermVelocityConvection_dynamic(gridGhost, tNow, SET, LLSid, cof, vyPrime_no_stencil)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
ns              = gridGhost.stencil;
%
row_no_stencil  = LLSid.rowLLS - ns;
%%
N               = length(LLSid.rowLLS);
vy              = zeros(1, N);
vx              = zeros(1, N);
%%%%
N_no_stencil    = gridGhost.N(1) - 2*ns;

% AO: The following was really slow. Replaced by an equivalent statement
% which is faster?
%for ss = 1 : N_no_stencil
% vy(row_no_stencil == ss) = vyPrime_no_stencil(ss);
%end
vy = vyPrime_no_stencil(row_no_stencil);

vy(row_no_stencil <= 0) = vyPrime_no_stencil(1);  % let the stencil value equal that at the base
vy(row_no_stencil > N_no_stencil) = vyPrime_no_stencil(N_no_stencil);  %
vy                  = vy + SET.vMean;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% cut off function: cof
Velocity{1}         = cof.*vy;
Velocity{2}         = cof.*vx;
%
end
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  end  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
function ydot = FcnTermCV_sL(derivL, derivR, V, sL, curvature, gradMag, SET)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% This function is used to calculate the time derivative of G
%
%%%%%%%%%%%%%%%%%%%%%%% Convection + sL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% D_t G = -(V(x,t) - sL \grad G/ \| \grad G \| ) \cdot \grad G - b k \|
% \grad G \|)

% b = -sL0*L
%
% Input: [derivL, derivR] are computed by the WENO method
% V is current flow velocity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim     = 2;
deriv   = cell(dim, 1);
for kk  = 1 : dim
    flowL = ((sL .* derivR{kk} <= 0) & (sL .* derivL{kk} <= 0));
    flowR = ((sL .* derivR{kk} >= 0) & (sL .* derivL{kk} >= 0));
    flows = ((sL .* derivR{kk} <  0) & (sL .* derivL{kk} >  0));
    if(any(flows(:)))
        conv = find(flows);
        s = zeros(size(flows));
        s(conv) = sL(conv) .* (abs(derivR{kk}(conv)) - abs(derivL{kk}(conv))) ...
            ./ (derivR{kk}(conv) - derivL{kk}(conv));
        
        flowL(conv) = flowL(conv) | (s(conv) < 0);
        flowR(conv) = flowR(conv) | (s(conv) >= 0);
    end
    deriv{kk} = derivL{kk} .* flowR + derivR{kk} .* flowL;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute magnitude of gradient.
mag   = zeros(size(derivL{1}));
for kk = 1 : dim;
    mag = mag + deriv{kk}.^2;
end
% mag represent \|\grad G\|
mag = max(sqrt(mag), eps);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vcomb   = cell(dim, 1);
for kk  = 1 : dim
    % Effective velocity field (for timestep bounding).
    Vcomb{kk} = V{kk} - sL .* deriv{kk} ./ mag;
end
dim     = 2;
delta   = zeros(size(derivL{1}));
% Determine the upwind direction dimension by dimension
for kk  = 1 : dim
    % Figure out upwind direction.
    flowL   = (Vcomb{kk} < 0);
    flowR   = (Vcomb{kk} > 0);
    deriv   = derivL{kk} .* flowR + derivR{kk} .* flowL;
    delta   = delta + deriv .* Vcomb{kk};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% item related to the curvature
% methane-air flame, T0 = 300 K, p0 = 1 bar, phi = 1.0
% flame speed
sL0     = SET.sL0;
% Markstein length
%     L       = 1e-3;
L       = 0; %AO
b       = -1*sL0*L;
delta   = delta + b .* curvature .* gradMag;
%
ydot    = -delta;
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataGhost = FcnLLS2dataGhost(dataLLS, dataGhost, LLSid, gridGhost)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% return the LLS data to the dataGhost data
%
%idxLLSGhost     = sub2ind(size(dataGhost), LLSid.rowLLS, LLSid.colLLS);
idxLLSGhost     = size(dataGhost,1)*(LLSid.colLLS-1) + LLSid.rowLLS;
dataGhost(idxLLSGhost) = dataLLS;

% The periodic ghost data need to be updated
%
indexPBaxis     = gridGhost.N(2) - gridGhost.stencil;  % index of the axis of the periodic boundary
dataGhost(:, indexPBaxis + 1 :  1 : indexPBaxis + gridGhost.stencil) = ...
    dataGhost(:, indexPBaxis - 1 : -1 : indexPBaxis - gridGhost.stencil);
%
% Flame base reattach
dataGhost(gridGhost.indexBaseY, gridGhost.indexBaseX) = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  breaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dataGhost = FcnLLS_ReINI_Global(dataGhost, gridGhost, NReINIGlobal, dtReINIGlobal)
%%
% global reinitialisation
% author: Jingxuan Li
% last updated: 2017-02-22
tReINI          = 0;
dataGhost0      = dataGhost;
for ss      = 1 : NReINIGlobal
    tSpan   = [tReINI, tReINI + dtReINIGlobal];
    [tReINI, dataGhost]   = FcnTVDRunge_Kutta3_ReINI_Global(tSpan, dataGhost, dataGhost0, gridGhost);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, y] = FcnTVDRunge_Kutta3_ReINI_Global(tSpan, y0, yRaw, gridGhost)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% y0 is the local level set data per substep
% yRaw is the local level set data before the reinitilisation step
t       = tSpan(1);
y       = y0;
ydot    = 0*y;
deltaT  = diff(tSpan);
Nst     = gridGhost.stencil;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydot(Nst + 1 : gridGhost.N(1) - Nst, Nst + 1 : gridGhost.N(2) - Nst)...
    = FcnLLS_ReINI_Global_ydot(y, yRaw, gridGhost);
%
t1      = t + deltaT;
y1      = y + deltaT*ydot;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydot(Nst + 1 : gridGhost.N(1) - Nst, Nst + 1 : gridGhost.N(2) - Nst)...
    = FcnLLS_ReINI_Global_ydot(y1, yRaw, gridGhost);
t2      = t1 + deltaT;
y2      = y1 + deltaT*ydot;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tHalf   = 0.25*(3*t + t2);
yHalf   = 0.25*(3*y + y2);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydot(Nst + 1 : gridGhost.N(1) - Nst, Nst + 1 : gridGhost.N(2) - Nst)...
    = FcnLLS_ReINI_Global_ydot(yHalf, yRaw, gridGhost);
tThreeHalf  = tHalf + deltaT;
yThreeHalf  = yHalf + deltaT*ydot;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%.
t       = (1/3)*(t + 2*tThreeHalf);
y       = (1/3)*(y + 2*yThreeHalf);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ydot = FcnLLS_ReINI_Global_ydot(dataGhost, dataGhost0, gridGhost)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% calculate the spatial derivative using WENO5 approach
[derivL,derivR] = FcnUpWindFirstWENO5_Global(dataGhost, gridGhost);
stencil = gridGhost.stencil;
G0      = dataGhost0(   stencil + 1 : gridGhost.N(1) - stencil,...
    stencil + 1 : gridGhost.N(2) - stencil);

ydot = FcnTermReINI_Global(derivL, derivR, G0, gridGhost);
end
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ydot = FcnTermReINI_Global(derivL, derivR, G0, gridGhost)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
%%%%%%%%%%%%%%%%%%%%%%% Normal term %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  D_t G =  S(x, t) (1 - \| \grad G \|)
%
% Input: [derivL, derivR] are computed by the WENO method
% sL is the current speed of sound, it may be constant or depend on time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FlagMethod = 3;
switch FlagMethod
    case 3
        dim = 2;
        %---------------------------------------------------------------------------
        % Compute Godunov derivative approximation for each dimension.
        deriv = cell(dim, 1);
        S     = G0;
        for kk = 1 : dim
            flowL = ((S .* derivR{kk} <= 0) & (S .* derivL{kk} <= 0));
            flowR = ((S .* derivR{kk} >= 0) & (S .* derivL{kk} >= 0));
            flows = ((S .* derivR{kk} <  0) & (S .* derivL{kk} >  0));
            if(any(flows(:)))
                conv = find(flows);
                s = zeros(size(flows));
                s(conv) = S(conv) .* (abs(derivR{kk}(conv)) - abs(derivL{kk}(conv))) ...
                    ./ (derivR{kk}(conv) - derivL{kk}(conv));
                flowL(conv) = flowL(conv) | (s(conv) < 0);
                flowR(conv) = flowR(conv) | (s(conv) >= 0);
            end
            deriv{kk} = derivL{kk} .* flowR + derivR{kk} .* flowL;
        end
        %---------------------------------------------------------------------------
        % Compute magnitude of gradient.
        mag   = zeros(size(derivL{1}));
        for kk = 1 : dim;
            mag = mag + deriv{kk}.^2;
        end
        mag = max(sqrt(mag), eps);
        %---------------------------------------------------------------------------
        dx          = min(gridGhost.dx);
        S           = G0./(G0.^2 + mag.^2.*dx.^2).^0.5;
        delta = -S;
        % Compute change in function and bound on step size.
        for kk = 1 : dim
            % Effective velocity field (for timestep bounding).
            v = S .* deriv{kk} ./ mag;
            delta = delta + v .* deriv{kk};
        end
end
ydot = -delta;
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [derivL,derivR] = FcnUpWindFirstWENO5_Global(dataGhost, gridGhost)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
%---------------------------------------------------------------------------
dxInv   = 1./gridGhost.dx;
dim     = 2;
%
D1     = cell(dim, 1);
derivL  = cell(dim, 1);
derivR  = cell(dim, 1);
stencil = gridGhost.stencil;
%
for kk = 1 : dim
    D1{kk}      = diff(dataGhost, 1, kk).*dxInv(kk);
    derivL{kk}  = FcnUpWindFirstWENO5_Global_One_Direction(D1{kk}, gridGhost.N, kk, 0, stencil);
    derivR{kk}  = FcnUpWindFirstWENO5_Global_One_Direction(D1{kk}, gridGhost.N, kk, 1, stencil);
end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function deriv = FcnUpWindFirstWENO5_Global_One_Direction(D1, N, idxDim, isRight, stencil)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% We always carry calculation in the column direction, if idxDim == 2, the
% matrix should be transposed
% N is the size of matrix not transposed
N0 = stencil - 3;
if (isRight == 1) && (idxDim == 1)
    D1 = D1(end:-1:1, :);
elseif (isRight == 1) && (idxDim == 2)
    D1 = D1(:, end:-1:1);
end
id1     = cell(5, 1);
id2     = cell(1, 1);
if  idxDim == 1
    % column direction
    temp1   = N0 + (1 : N(1) - 2*stencil);
    id2{1}  = stencil + 1 : N(2) - stencil;
else
    % row direction
    temp1   = N0 + (1 : N(2) - 2*stencil);
    id2{1}  = stencil + 1 : N(1) - stencil;
    D1      = D1';
end
for kk  = 1 : 5
    id1{kk} = temp1 + kk - 1;
end
% The third order accurate ENO scheme
phi{1}  =  1/3*D1(id1{1}, id2{1}) - 7/6*D1(id1{2}, id2{1}) + 11/6*D1(id1{3}, id2{1});
phi{2}  = -1/6*D1(id1{2}, id2{1}) + 5/6*D1(id1{3}, id2{1}) +  1/3*D1(id1{4}, id2{1});
phi{3}  =  1/3*D1(id1{3}, id2{1}) + 5/6*D1(id1{4}, id2{1}) -  1/6*D1(id1{5}, id2{1});
%---------------------------------------------------------------------------
% Smoothness estimates.
smooth = cell(3,1);
smooth{1} = (13/12)*(  D1(id1{1}, id2{1}) - 2*D1(id1{2}, id2{1}) +   D1(id1{3}, id2{1})).^2 ...
    +  (1/4)*(  D1(id1{1}, id2{1}) - 4*D1(id1{2}, id2{1}) + 3*D1(id1{3}, id2{1})).^2;
smooth{2} = (13/12)*(  D1(id1{2}, id2{1}) - 2*D1(id1{3}, id2{1}) +   D1(id1{4}, id2{1})).^2 ...
    +  (1/4)*(  D1(id1{2}, id2{1})                        -   D1(id1{4}, id2{1})).^2;
smooth{3} = (13/12)*(  D1(id1{3}, id2{1}) - 2*D1(id1{4}, id2{1}) +   D1(id1{5}, id2{1})).^2 ...
    +  (1/4)*(3*D1(id1{3}, id2{1}) - 4*D1(id1{4}, id2{1}) +   D1(id1{5}, id2{1})).^2;
%
%--------------------------------------------------------------------------
% weight
weight = [0.1; 0.6; 0.3];
%
%--------------------------------------------------------------------------
% epsilon
epsilonCalMethod = 'maxOverGrid';
switch(epsilonCalMethod)
    case 'constant'
        epsilon     = 1e-6;
    case 'maxOverGrid'
        D1squared   = D1.^2;
        epsilon     = 1e-6*max(D1squared(:)) + 1e-99;
end
%
%---------------------------------------------------------------------------
% Compute and apply weights to generate a higher order WENO approximation.
deriv = weightWENO(phi, smooth, weight, epsilon);
if  idxDim == 2
    deriv = deriv';
end
if (isRight == 1) && (idxDim == 1)
    deriv = deriv(end:-1:1, :);
elseif (isRight == 1) && (idxDim == 2)
    deriv = deriv(:, end:-1:1);
end
end
%
%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function deriv = weightWENO(d, s, w, epsilon)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% Compute weighting terms
alpha1  = w(1)./(s{1} + epsilon).^2;
alpha2  = w(2)./(s{2} + epsilon).^2;
alpha3  = w(3)./(s{3} + epsilon).^2;
sum     = (alpha1 + alpha2 + alpha3);
deriv   = (alpha1.*d{1} + alpha2.*d{2} + alpha3.*d{3})./sum;
end
%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  breaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dataLLS = FcnLLS_ReINI(dataLLS0, dataGhost, gridGhost, LLSid, NReINI, dtReINI)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% reinitialisation
tReINI      = 0;
dataLLS     = dataLLS0;
dataGhost   = FcnLLS2dataGhost(dataLLS, dataGhost, LLSid, gridGhost);
clear FF
FF = 0;
% dtReINI = dt;
for ss      = 1 : NReINI
    tSpan   = [tReINI, tReINI + dtReINI];
    [tReINI, dataLLS]   = FcnTVDRunge_Kutta3_ReINI(tSpan, dataLLS, dataLLS0, dataGhost, gridGhost, LLSid, FF);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, y] = FcnTVDRunge_Kutta3_ReINI(tSpan, y0, yRaw, dataGhost, gridGhost, LLSid, FF)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
t       = tSpan(1);
y       = y0;
deltaT  = diff(tSpan);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydot    = FcnLLS_ReINI_ydot(y, yRaw, dataGhost, gridGhost, LLSid, FF);
%
t1      = t + deltaT;
y1      = y + deltaT*ydot;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydot    = FcnLLS_ReINI_ydot(y1, yRaw, dataGhost, gridGhost, LLSid, FF);
t2      = t1 + deltaT;
y2      = y1 + deltaT*ydot;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tHalf   = 0.25*(3*t + t2);
yHalf   = 0.25*(3*y + y2);
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ydot    = FcnLLS_ReINI_ydot(yHalf, yRaw, dataGhost, gridGhost, LLSid, FF);
tThreeHalf  = tHalf + deltaT;
yThreeHalf  = yHalf + deltaT*ydot;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t       = (1/3)*(t + 2*tThreeHalf);
y       = (1/3)*(y + 2*yThreeHalf);
%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ydot = FcnLLS_ReINI_ydot(dataLLS, dataLLS0, dataGhost, gridGhost, LLSid, FF)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% dataLLS is the local level set data per substep
% calculate the spatial derivative using WENO5 approach
% update the dataGhost
dataGhost       = FcnLLS2dataGhost(dataLLS, dataGhost, LLSid, gridGhost);
% construct the matrix for the spatial derivative calculation
Mtr             = FcnLLS_MtrWENO(dataGhost, LLSid);
% calculate the spatial derivative
[derivL,derivR] = FcnUpWindFirstWENO5(Mtr, gridGhost);
% calculate the forcing term
ydot            = FcnTermReINI_Peng(derivL, derivR, dataLLS0, gridGhost, 0);
end
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
function [derivL,derivR] = FcnUpWindFirstWENO5(Mtr, gridGhost)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
derivL = cell(2,1);
derivR = cell(2,1);
for dim = 1 : 2
    dxInv = 1/gridGhost.dx(dim);
    %
    Array   = Mtr{dim};
    D1L     = dxInv*diff(Array, 1);
%     D1R     = flipud(D1L);
    D1R     = D1L(end:-1:1,:);
    %
    derivL{dim} = FcnUpWindFirstWENO5_one_direction(D1L);
    derivR{dim} = FcnUpWindFirstWENO5_one_direction(D1R);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% start %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function deriv = FcnUpWindFirstWENO5_one_direction(D1)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% The third order accurate ENO scheme
phi{1} =  1/3*D1(1,:) - 7/6*D1(2,:) + 11/6*D1(3,:);
phi{2} = -1/6*D1(2,:) + 5/6*D1(3,:) +  1/3*D1(4,:);
phi{3} =  1/3*D1(3,:) + 5/6*D1(4,:) -  1/6*D1(5,:);
%---------------------------------------------------------------------------
% Smoothness estimates.
smooth = cell(3,1);
smooth{1} = (13/12)*(  D1(1,:) - 2*D1(2,:) +   D1(3,:)).^2 ...
    +  (1/4)*(  D1(1,:) - 4*D1(2,:) + 3*D1(3,:)).^2;
smooth{2} = (13/12)*(  D1(2,:) - 2*D1(3,:) +   D1(4,:)).^2 ...
    +  (1/4)*(  D1(2,:) -   D1(4,:)).^2;
smooth{3} = (13/12)*(  D1(3,:) - 2*D1(4,:) +   D1(5,:)).^2 ...
    +  (1/4)*(3*D1(3,:) - 4*D1(4,:) +   D1(5,:)).^2;
%
%--------------------------------------------------------------------------
% weight
weight = [0.1; 0.6; 0.3];
%
%--------------------------------------------------------------------------
% epsilon
epsilonCalMethod = 'maxOverGrid';
switch(epsilonCalMethod)
    case 'constant'
        epsilon     = 1e-6;
    case 'maxOverGrid'
        D1squared   = D1.^2;
        epsilon     = 1e-6*max(D1squared(:)) + 1e-99;
end
%
%---------------------------------------------------------------------------
% Compute and apply weights to generate a higher order WENO approximation.
deriv = weightWENO(phi, smooth, weight, epsilon);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function ydot = FcnTermReINI_Peng(derivL, derivR, G0, gridGhost, F)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
%
% We use the smoothed sign function S to identify the flow direction, this
% term is considered to be normal to the flame front and is a scalar.
%
%%%%%%%%%%%%%%%%%%%%%%% Normal term %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  D_t G =  S(x, t) (1 - \| \grad G \|) + beta F
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = 2;
%---------------------------------------------------------------------------
% Compute Godunov derivative approximation for each dimension.
deriv = cell(dim, 1);
S     = G0;
for kk = 1 : dim
    flowL = ((S .* derivR{kk} <= 0) & (S .* derivL{kk} <= 0));
    flowR = ((S .* derivR{kk} >= 0) & (S .* derivL{kk} >= 0));
    flows = ((S .* derivR{kk} <  0) & (S .* derivL{kk} >  0));
    if(any(flows(:)))
        conv = find(flows);
        s = zeros(size(flows));
        s(conv) = S(conv) .* (abs(derivR{kk}(conv)) - abs(derivL{kk}(conv))) ...
            ./ (derivR{kk}(conv) - derivL{kk}(conv));
        flowL(conv) = flowL(conv) | (s(conv) < 0);
        flowR(conv) = flowR(conv) | (s(conv) >= 0);
    end
    deriv{kk} = derivL{kk} .* flowR + derivR{kk} .* flowL;
end
%---------------------------------------------------------------------------
% Compute magnitude of gradient.
mag   = zeros(size(derivL{1}));
for kk = 1 : dim;
    mag = mag + deriv{kk}.^2;
end
mag = max(sqrt(mag), eps);
%---------------------------------------------------------------------------
dx          = min(gridGhost.dx);
S           = G0./(G0.^2 + mag.^2.*dx.^2).^0.5;
delta = -S;
for kk = 1 : dim
    v = S .* deriv{kk} ./ mag;
    delta = delta + v .* deriv{kk};
end
ydot = -delta;
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  breaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mtr = FcnLLS_MtrWENO(dataGhost, LLSid)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% This function is used to construct the matrix used in the WENO approach
dataX           = zeros(7, length(LLSid.rowLLS));
dataY           = zeros(7, length(LLSid.rowLLS));
% rowLLS and colLLS now correspond to the indices of the ghost data
rowLLSGhost     = LLSid.rowLLS;
colLLSGhost     = LLSid.colLLS;
sz              = size(dataGhost,1);
for kk = 1 : 7
%     idxY        = sub2ind(size(dataGhost), rowLLSGhost + kk - 4, colLLSGhost);
%     idxX        = sub2ind(size(dataGhost), rowLLSGhost, colLLSGhost + kk - 4);
    idxY        = rowLLSGhost + kk - 4 + (colLLSGhost          -1)*sz;
    idxX        = rowLLSGhost          + (colLLSGhost + kk - 4 -1)*sz;
    dataX(kk,:) = dataGhost(idxX);
    dataY(kk,:) = dataGhost(idxY);
end
Mtr{1}  = dataY;
Mtr{2}  = dataX;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%  breaks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [LLSid, dataLLS] = FcnLLS_INI(dataGhost, gridGhost)
%%
% author: Jingxuan Li
% last updated: 2017-02-22
% The method proposed in Peng et al. (1999) is used to identify the
% different zones used in the local level set method.
stencil = gridGhost.stencil;
% Number of cells to describe the halo zone + tube zone, tube zone
data    = dataGhost((1 + stencil : gridGhost.N(1) - stencil), (1 + stencil : gridGhost.N(2) - stencil));
data_abs = abs(data);
%%%%%%%%%%%%%%%%%%%%%%%%%% start of the code %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
gamma   = [9, 6, 3].*max( gridGhost.dx );
%
[rowM4, colM4]  = find(data_abs <  gamma(3));
[rowM3, colM3]  = find(data_abs >= gamma(3) & data_abs < gamma(2));
%AO: Change -1
%[rowM2, colM2]  = find(data_abs >= gamma(2) & data_abs < gamma(1));
rowM2 = [];
colM2 = [];
%B       = data_abs < gamma(1);
% AO Change 0
B       = data_abs < gamma(2);
C1      = diff(B,1,1);
C2      = diff(B,1,2);
NData   = gridGhost.N - 2*stencil;
BdiffyM = cat(1, C1, zeros(1, NData(2)));
BdiffyP = cat(1, zeros(1, NData(2)), C1);
BdiffxM = cat(2, zeros(NData(1), 1), C2);
BdiffxP = cat(2, C2, zeros(NData(1), 1));
E       = (B == 0 & (BdiffxP ~= 0 | BdiffxM ~=0 | BdiffyP ~= 0 | BdiffyM ~= 0));
[rowM1, colM1]  = find(E ~= 0);
NMargin(1)  = length(rowM2) + length(rowM3) + length(rowM4);
NMargin(2)  = length(rowM3) + length(rowM4);
NMargin(3)  = length(rowM4);
rowLLS  = cat(1, rowM4, rowM3, rowM2, rowM1) + stencil;
colLLS  = cat(1, colM4, colM3, colM2, colM1) + stencil;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%idxLLS      = sub2ind(size(dataGhost), rowLLS, colLLS);
idxLLS      = size(dataGhost,1)*(colLLS-1) + rowLLS;
dataLLS     = dataGhost(idxLLS);
sizeDataLLS = size(dataLLS);
if (sizeDataLLS(1) >= sizeDataLLS(2))
    dataLLS  = dataLLS';
end
LLSid.rowLLS    = rowLLS;
LLSid.colLLS    = colLLS;
LLSid.NMargin   = NMargin;
LLSid.gamma     = gamma;
%
end
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


