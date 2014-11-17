function uRatio = Fcn_TD_calculation_uRatio(Var,indexSection)
% This function is used to calculate the velocity ratio before the flame,
% uRatio = (AP(t-tauP-tauf) - AM(t-tauf))/(rho*c*u);
% The time delay of the flame describing function is included.
% index of time steps Var(1):Var(2)
% indexSection denotes the index of section before the flame
% last edited: 2011-11-12 16:04
global CI
ss  = indexSection;
AP  = Fcn_interp1_varied_td(    CI.TD.AP(ss,:),...
                                Var,...
                                CI.TD.tauf(Var(1):Var(2)) + CI.TP.tau_plus(ss),...
                                CI.TD.dt);
AM  = Fcn_interp1_varied_td(    CI.TD.AM(ss,:),...
                                Var,...
                                CI.TD.tauf(Var(1):Var(2)),...
                                CI.TD.dt);
uRatio = (AP - AM)./CI.TP.RhoCU(ss);
% CI.TP.RhoCU denotes rho*c*u
% velocity ratio limit : CI.TD.uRatioMax, which is first defined in the
% program : Fcn_TD_initialization
uRatio = Fcn_saturation(uRatio,CI.TD.uRatioMax);
%
% ----------------------------end------------------------------------------