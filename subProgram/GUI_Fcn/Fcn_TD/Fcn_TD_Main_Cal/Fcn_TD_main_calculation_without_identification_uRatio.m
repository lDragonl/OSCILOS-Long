function Fcn_TD_main_calculation_without_identification_uRatio
global CI
%
% --------------------------------------
hWaitBar = waitbar(0,'Time domain calculations, please wait...');
%
%
for nn = 1:CI.TD.nRound
    Var(1:2)     = CI.TD.nPadding + [(nn-1)*CI.TD.nGap + 1,...
                                 nn*CI.TD.nGap];                     
%      if CI.FM.NL.style == 4 % G-equation case
%          CI.FM.NL.Model4.IT = nn;
%      end
    Fcn_TD_calculation_one_gap(Var)
    waitbar(nn/CI.TD.nRound);
    drawnow
end
close(hWaitBar)
assignin('base','CI',CI);
% Fcn_TD_results_plot                                     % Plot the result
%
% ----------------------------end------------------------------------------