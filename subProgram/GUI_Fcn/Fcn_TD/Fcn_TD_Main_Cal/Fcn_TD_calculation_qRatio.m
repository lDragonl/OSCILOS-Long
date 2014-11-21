function qRatio = Fcn_TD_calculation_qRatio(uRatio,yFTFGreen,Var,dt,Lf)
% This function is used to calculate qRatio of flame
% uRatio is the velocity ratio before the flame
% yFTFGreen is the Green's function of the flame transfer function for weak
% flow perturbations
% index of time steps Var(1):Var(2)
% dt is the time step
% Lf is the nonlinear gain ratio
% time delay of flame transfer function is inclued in uRatio
%
% last edited: 2014-11-12 16:11
%
global CI
NGreen  = length(yFTFGreen);
index   = Var(1)-NGreen+1 : Var(2);
qRatio  = conv(uRatio(index),yFTFGreen,'valid').*dt;
switch CI.FM.NL.style
    case 2
        qRatio = Fcn_saturation(qRatio,CI.FM.NL.Model2.alpha);
    case 4 % Case of the G-equation
        error(' not implemented yet')
    otherwise
        qRatio = qRatio.*Lf((Var(1):Var(2)));
%         qRatio = qRatio.*0.2;
end
%
% --------------------------------end--------------------------------------
