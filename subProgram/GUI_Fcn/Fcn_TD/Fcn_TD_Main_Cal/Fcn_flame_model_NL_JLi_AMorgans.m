function [Lf,tauf] = Fcn_flame_model_NL_JLi_AMorgans(uRatio)
global CI
% This function is used to calculate the nonlinear gain ratio and time
% delay of flame describing function of our model
%
% last edit: 2014-11-13 
%
Lf      = interp1(  CI.FM.NL.Model3.uRatio,...
                    CI.FM.NL.Model3.Lf,...
                    uRatio,'linear','extrap');
tauf    = CI.FM.FTF.tauf + CI.FM.NL.Model3.taufN.*(1-Lf);
%
% ---------------------------end-------------------------------------------
