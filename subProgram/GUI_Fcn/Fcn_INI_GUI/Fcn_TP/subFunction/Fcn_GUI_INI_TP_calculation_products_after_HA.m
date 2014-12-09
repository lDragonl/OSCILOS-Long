function [TRatio, cMean2, DeltaHr, gamma2, Cp2] =...
         Fcn_GUI_INI_TP_calculation_products_after_HA(handles,indexHA_num,TMean1,pMean1)
% This function is used to calculate the mean properties after the heat
% addition interface
%
% indexHA_num is the order of heat addition interface from the inlet to the
% outlet 
% TMean1 denotes the incident mean temperature, 
% pMean1 denotes the incident mean pressure
% TRatio is the temperature jump ratio
% cMean2 denotes the speed of sound after the interface
% DeltaHr denotes the heat addtion
% gamma2 denotes the specific heat ratio
% Cp2 denotes the heat capacity at constant pressure after the interface
%
% first created: 2014-12-03
% last modified: 2014-12-03
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
global CI
ss          = indexHA_num;
HA_style    = CI.TP.HA_style(ss);                           % heat addition style of selected heating interface
switch HA_style
    case 1
        % temperature determined
        TRatio  = CI.TP.TRatio(ss);
        TMean2  = TMean1*TRatio;
        [temp1,cMean2,DeltaHr,Cp2] = Fcn_calculation_c_q_air(TMean1,TMean2);
        gamma2  = Cp2./(Cp2-CI.R_air);     
    case 2
        % fuel determined
%         CI.TP.indexFuel = get(handles.pop_FD_fuel,'value');                     % index of fuel
%         CI.TP.eff       = str2num(get(handles.edit_FD_effi,     'string'));     % combustion efficiency
%         CI.TP.Phi       = str2num(get(handles.edit_FD_phi,      'string'));     % equivalence ratio
%         CI.TP.dil       = str2num(get(handles.edit_FD_dilute,   'string'));     % diluted ratio
    [   TMean2,...
        chi,...
        DeltaHr,...
        cMean2,...
        Cp2,...
        WProd ]   = ...
            Fcn_calculation_combustion_products(CI.TP.indexFuel(ss),...
                                                CI.TP.Phi(ss),...
                                                TMean1,...
                                                pMean1,...
                                                CI.TP.eff(ss),...
                                                CI.TP.dil(ss));
        TRatio              = TMean2./TMean1;                           % Temperature jump ratio
        RProd               = CI.Ru./(WProd).*1000;                     % gas constant per mass for the combustion product
        gamma2              = Cp2./(Cp2-RProd);
    otherwise
        % Code for when there is no match.
end
assignin('base','CI',CI)
% -------------------------------------------------------------------------