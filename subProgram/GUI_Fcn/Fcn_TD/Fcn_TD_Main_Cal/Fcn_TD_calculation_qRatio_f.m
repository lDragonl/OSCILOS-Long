function qRatio = Fcn_TD_calculation_qRatio_f(uRatio,yFTFGreen,Var,dt,Lf)
% The flame transfer function is the function of f, but s. 
% This function is used to calculate qRatio of flame
% uRatio is the velocity ratio before the flame
% yFTFGreen is the Green's function of the flame transfer function for weak
% flow perturbations
% index of time steps Var(1):Var(2)
% dt is the time step
% Lf is the nonlinear gain ratio
% time delay of flame transfer function is inclued in uRatio
%
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
% first created:    2014-11-11
% last edited:      2014-11-19
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
    case 5 % AO: Case of Convective G-equation
        error('This can be done but needs to be coded. See comments for more info');
        % This can be done in two ways:
        %1) There is an analytical expression for the FTF. I can code it
        %easily
        %2) One can even run time domain simulations to calculate the FTF
        %(useless because there is the analytical model, I have tested this
        %already and the two methods work) OR to calculate the FDF. This is
        %more useful probably. Need to check with Aimee.
    otherwise
        qRatio = qRatio.*Lf((Var(1):Var(2)));
%         qRatio = qRatio.*0.2;
end
%
% --------------------------------end--------------------------------------
