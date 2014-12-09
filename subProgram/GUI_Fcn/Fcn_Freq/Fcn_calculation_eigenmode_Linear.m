function [x_resample,p,u]=Fcn_calculation_eigenmode_Linear(s_star)
% This function is used to plot the modeshape of slected mode, 
% Linear
global CI
[R1,R2]     = Fcn_boundary_condition(s_star);
Rs          = -0.5*CI.TP.M_mean(end)./(1 + 0.5*(CI.TP.gamma(end) - 1 ).*CI.TP.M_mean(end));
%
%
Te          = Fcn_TF_entropy_convection(s_star);
%--------------------------------
tau_plus    = CI.TP.tau_plus;
tau_minus   = CI.TP.tau_minus;
tau_c       = CI.TP.tau_c;
%
% -------------------------------------------------------------------------
A_minus(1)  = 1;
E(1)        = 0;
A_plus(1)   = R1.*A_minus(1);
Array(:,1)  = [A_plus(1),A_minus(1),E(1)]';
%

% -------------------------------------------------------------------------
for ss = 1:CI.TP.numSection-1 
    D1 = diag([ exp(-s_star*tau_plus(ss)),...
                exp( s_star*tau_minus(ss)),...
                exp(-s_star*tau_c(ss))]);
    switch CI.CD.SectionIndex(ss+1)
        case 0
            CI.TPM.Z{ss}       = CI.TPM.BC{ss}*D1;
        case 10
            B2b     = zeros(3);
            B2b(3,2)= 0;
            Bsum    = CI.TPM.B1{2,ss}*(CI.TPM.B2{1,ss}\CI.TPM.B1{1,ss}) + B2b;
            BC1     = Bsum*CI.TPM.C1;
            BC2     = CI.TPM.B2{2,ss}*CI.TPM.C2;
            CI.TPM.Z{ss}       = (BC2\BC1)*D1;
        case 11
            % linear flame transfer function
            FTF     = Fcn_flame_model(s_star);
            B2b     = zeros(3);
            temp    = D1*Array(:,ss);
            uf      = abs((temp(1) - temp(2))./(CI.TP.c_mean(1,ss).*CI.TP.rho_mean(1,ss)));         % velocity before the flame
            B2b(3,2)= CI.TP.DeltaHr./CI.TP.c_mean(2,ss+1)./CI.TP.c_mean(1,ss)./CI.TP.Theta(ss).*FTF;
            Bsum    = CI.TPM.B1{2,ss}*(CI.TPM.B2{1,ss}\CI.TPM.B1{1,ss}) + B2b;
            BC1     = Bsum*CI.TPM.C1;
            BC2     = CI.TPM.B2{2,ss}*CI.TPM.C2;
            CI.TPM.Z{ss}       = (BC2\BC1)*D1;
    end
    Array(:,ss+1)    = CI.TPM.Z{ss}*Array(:,ss);
end
%
for k = 1:CI.TP.numSection-1
    A_plus(k+1)     = Array(1,k+1);
    A_minus(k+1)    = Array(2,k+1);
    E(k+1)          = Array(3,k+1);
end
%
for k=1:length(CI.CD.x_sample)-1
    x_resample(k,:)=linspace(CI.CD.x_sample(k),CI.CD.x_sample(k+1),100);    % resample of x-coordinate
end
% speed of sound and mean velocity
c_mean      = CI.TP.c_mean(1,:);
u_mean      = CI.TP.u_mean(1,:);
rho_mean    = CI.TP.rho_mean(1,:);
%     
for k = 1:length(CI.CD.x_sample)-1
    kw1_plus(k) = s_star./(c_mean(k)+u_mean(k));
    kw1_minus(k) = s_star./(c_mean(k)-u_mean(k));
    p(k,:) =    A_plus(k).*exp(-kw1_plus(k).*(x_resample(k,:)-x_resample(k,1)))+...
                A_minus(k).*exp( kw1_minus(k).*(x_resample(k,:)-x_resample(k,1)));
    u(k,:) =    (A_plus(k).*exp(-kw1_plus(k).*(x_resample(k,:)-x_resample(k,1)))-...
                A_minus(k).*exp( kw1_minus(k).*(x_resample(k,:)-x_resample(k,1))))./rho_mean(k)./c_mean(k);
end


%
%----------------------Pressure Reflection coefficients -------------------
%
function [R1,R2] = Fcn_boundary_condition(s)
global CI
R1      = polyval(CI.BC.num1,s)./polyval(CI.BC.den1,s).*exp(-CI.BC.tau_d1.*s);
R2      = polyval(CI.BC.num2,s)./polyval(CI.BC.den2,s).*exp(-CI.BC.tau_d2.*s);
%
%----------------------Entropy convection transfer function ---------------
%
function Te = Fcn_TF_entropy_convection(s)
global CI
tau     = CI.BC.ET.Dispersion.Delta_tauCs;
k       = CI.BC.ET.Dissipation.k;
switch CI.BC.ET.pop_type_model
    case 1
        Te = 0;
    case 2
        Te = k.*exp((tau.*s).^2./4);
    case 3
        if tau == 0
            tau = eps;
        end
        Te = k.*(exp(tau*s) - exp(-tau*s))./(2*tau);
end                        
%
% ----------------------linear Flame transfer function --------------------
%
% !!! further development
% function F = Fcn_flame_model(s,indexHP)
% global CI
% num = CI.FM.FTF.num{indexHP};
% den = CI.FM.FTF.den{indexHP};
% tauf = CI.FM.FTF.tauf(indexHP);
% F = polyval(num,s)./polyval(den,s).*exp(-s.*tauf);            
function F = Fcn_flame_model(s)
global CI
num = CI.FM.FTF.num;
den = CI.FM.FTF.den;
tauf = CI.FM.FTF.tauf;
F = polyval(num,s)./polyval(den,s).*exp(-s.*tauf);        
%
% -----------------------------end-----------------------------------------