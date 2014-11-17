function Fcn_TD_main_calculation_period_uRatio_envelope(mm)
% calculate the envelope of a period
% last edit: 2014-11-13
%
global CI
Var(1)     = CI.TD.nPadding + (CI.TD.indexPeriod(mm,1)-1)*CI.TD.nGap + 1;
Var(2)     = CI.TD.nPadding +  CI.TD.indexPeriod(mm,2)*CI.TD.nGap;
N = Var(2) - Var(1) + 1;
CI.TD.uRatioEnv(Var(1):Var(2))...
    = Fcn_calculation_envelope(CI.TD.uRatio(Var(1):Var(2)),N);
% 
% [Lf,tauf] = Fcn_flame_model_NL_JLi_AMorgans(uRatio)
%
% ------------------------------end----------------------------------------
