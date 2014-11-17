function Fcn_TD_initialization_flame_model
% -------------------------------------------------------------------------
% ---------------------flame transfer function-----------------------------
% -------------------------------------------------------------------------
% For the first loop, the time delay and nonlinear gain ratio are set to
% the time delay of flame transfer function for weak perturbations and unit,
% respectively
%
% Last edit: 2014-11-13 16:22
%
global CI
uRatio0 = 0;                                                                % for weak velocity perturbations
switch  CI.FM.NL.style
    case 3
        [Lf,tauf] = Fcn_flame_model_NL_JLi_AMorgans(uRatio0);
    otherwise
        Lf = 1;
        tauf = CI.FM.FTF.tauf;
end
CI.TD.tauf = zeros(1,CI.TD.nTotal) + tauf;                                  % Nonlinear time delay, which varies with velocity perturbations
CI.TD.Lf   = zeros(1,CI.TD.nTotal) + Lf;                                    % Nonlinear model, which describes the saturation of heat release
                                                                            % rate with velocity perturbations
%
CI.TD.indexFlame = find(CI.CD.index == 1);
%
assignin('base','CI',CI);
%
% --------------------------end--------------------------------------------