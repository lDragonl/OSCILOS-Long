function Fcn_TD_main_calculation
global CI
%
% In case the flame describing function is determined by the envelope of
% velocity ratio, uRatio, iteration-convergence is necessary. 
% -------------------------------------------------------------------------
hWaitBar = waitbar(0,'Time domain calculations, please wait...');          
%
for mm = 1:CI.TD.nPeriod
    waitbar(mm/CI.TD.nPeriod);
    drawnow
    switch CI.FM.NL.style    
        case 3
            Iteration.tol = 1e-5;
            Iteration.num = 50;
            Fcn_TD_main_calculation_iteration_convergence(mm,Iteration)      
        case 4 % G-equation
            % Compute the flame shape and area
            [CI.FM.NL.Model4.q_ratio,CI.FM.NL.Model4.xi] = Fcn_TD_Gequ_interface...
                ( CI.FM.NL.Model4.SU, CI.FM.NL.Model4.xi, CI.FM.NL.Model4.y_vec, CI.TD.dt, 0, CI.FM.NL.Model4.U1, ...
                CI.FM.NL.Model4.area_ratio, uRatioEnv1,CI.TP.Q * CI.FM.NL.Model4.area_ratio, CI.TP.DeltaHr,CI.FM.NL.Model4.rho1 );
            
            % --------------------------
            Fcn_TD_calculation_one_period(mm)
            % --------------------------
            Fcn_TD_main_calculation_period_uRatio_envelope(mm)   % calculate the envelope                            
        otherwise
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



