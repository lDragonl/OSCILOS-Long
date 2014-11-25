function Fcn_TD_main_calculation
global CI
%
% In case the flame describing function is determined by the envelope of
% velocity ratio, uRatio, iteration-convergence is necessary. 
% -------------------------------------------------------------------------
hWaitBar = waitbar(0,'Time domain calculations, please wait...');
%
if CI.FM.NL.style ==4 % Convergence method notr required for G-EQuation
    
else
    for mm = 1:CI.TD.nPeriod
        switch CI.FM.NL.style
            case {3}
                Fcn_TD_main_calculation_iteration_convergence(mm)
                
            otherwise
                waitbar(mm/CI.TD.nPeriod);
                drawnow
                % --------------------------
                Fcn_TD_calculation_one_period(mm)
                % --------------------------
                Fcn_TD_main_calculation_period_uRatio_envelope(mm)   % calculate the envelope
        end
    end
close(hWaitBar)
assignin('base','CI',CI);

% 
figure
plot(CI.TD.tSpTotal,CI.TD.AP(1,:))

figure
plot(CI.TD.tSpTotal,CI.TD.uRatio,'-k')
hold on
plot(CI.TD.tSpTotal,CI.TD.uRatioEnv,'-r')


figure
plot(CI.TD.Lf)
%



