function [tauMin,taufMin] = Fcn_TD_minimum_time_delay(uRatioMax)
% This function is used to calculate the minimum time delay 
% In this version, the convection of entropy is not accounted for
% the function ``Fcn_PreProcessing '' must be run before this program
% last edit: 2014-11-12
%
global CI
% --------------------------------
% Convection time of the acoustic waves in different sections
tauMin = min(   min(CI.TP.tau_plus),...
                min(CI.TP.tau_minus));
%
% --------------------------------
% Boundary condition
% Inlet boundary:
tauMin = min(tauMin, CI.TP.tau_minus(1) + CI.BC.tau_d1);
% Outlet boundary:
tauMin = min(tauMin, CI.TP.tau_plus(1)  + CI.BC.tau_d2);
%
% --------------------------------
% nonlinear flame describing function time delay
switch CI.FM.NL.style
    case 1
        taufMin = CI.FM.FTF.tauf;
    case 2
        taufMin = CI.FM.FTF.tauf;
    case 3
        % Flame describing function model proposed by J.LI
        % and A.Morgans JSV
        % We must define a maximum velocity ratio value
        uRatioLimit     = [0 uRatioMax];
        LfLimit         = interp1(  CI.FM.NL.Model3.uRatio,...
                                    CI.FM.NL.Model3.Lf,...
                                    uRatioLimit,'linear','extrap');
        taufNLimit      = CI.FM.NL.Model3.taufN.*(1 - LfLimit);
        taufLimit       = CI.FM.FTF.tauf + taufNLimit;
        taufMin         = min(taufLimit);
end
tauMin          = min(tauMin, taufMin);
%
% ------------------------ end --------------------------------------------