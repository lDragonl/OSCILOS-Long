function Fcn_TD_calculation_one_gap(Var)
% This function is used to calculate the wave values in one gap
% one gap means nGap number of time steps, which are calculated in one loop
% The calculation is from the left boundary to the right boundary
% convolutions take most calculation resources
%
% author: Jingxuan LI (jingxuan.li@imperial.ac.uk)
% first created:    2014-11-11
% last edited:      2014-11-19
%
% AP is the pressure wave propagating in teh direction of flow
% AM is the pressure wave propagating in the opposite direction of flow
% the lines of AP and AM indicate the section (1 is upstream boundary, N is
% downstream boundary). The columns are the values as a function fo time. 
%
global CI

% from left to right
% -------------------------------------------------------------------------
% --------------------- Left boundary -------------------------------------
% -------------------------------------------------------------------------
% A1+(t) = conv(R1,A1-(t-tau1-))
CI.TD.AP(1, Var(1):Var(2))...
    = Fcn_lsim_varied_td(   CI.TD.AM(1,:),...
                            CI.TD.Green.sys{1},...
                            CI.TD.Green.tauConv2(1),...
                            Var,...
                            CI.TP.tau_minus(1) + CI.BC.tau_d1,...
                            CI.TD.dt);
CI.TD.AP(1, Var(1):Var(2)) = CI.TD.AP(1, Var(1):Var(2)) + CI.TD.pNoiseBG(Var(1):Var(2)); % add additional noise to wave propagating in direction of flow
%
% -------------------------------------------------------------------------
% ------------------ interfaces between two sections ----------------------
% -------------------------------------------------------------------------
for ss = 1:CI.TP.numSection-1 
    y(1,1:CI.TD.nGap)   = Fcn_interp1_varied_td(CI.TD.AP(ss,:),    Var, CI.TP.tau_plus(ss),    CI.TD.dt);
    y(2,1:CI.TD.nGap)   = Fcn_interp1_varied_td(CI.TD.AM(ss+1,:),  Var, CI.TP.tau_minus(ss+1), CI.TD.dt);
    y(3,1:CI.TD.nGap)  = 0; 
    switch CI.CD.index(ss+1)
        case 0   % only area change
            % ---------------------------
            %     [ A2+ ]       [ A1+ ]
            %     [ A1- ] =  [Z][ A2- ]
            %     [ E2  ]       [ E1  ]
            %
            x = CI.TD.IF.Z{ss}*y;
            % ---------------------------
        case 1   % with heat addition
            % ---------------------------
            %     [ A2+ ]       [ A1+ ]   [    ]
            %     [ A1- ] =  [Z][ A2- ] + [ Ar ] conv(FTF,u*rho*c)
            %     [ E2  ]       [ E1  ]   [    ]
            %     u*(rho*c) = A1+(t-tauf) - A1-(t-tauf)
            %
            if CI.FM.NL.style == 4 % In the case of the G-equation, the value of qratio is computed from a pde, not from a transfer function using Greens function
                
                CI.TD.uRatio(Var(1):Var(2)) = Fcn_TD_calculation_uRatio(Var,ss); % compute the values of uratio in every section
                
%                 amp = 0.17;
%                 freq = 333;
%                 alpha = 2;
%                 beta = 5.4;
%                 uratio = amp * (alpha * cos(freq* CI.TD.tSpTotal(Var(1):Var(2))) + beta * sin(freq * CI.TD.tSpTotal(Var(1):Var(2))));
%                 CI.TD.uRatio(Var(1):Var(2)) = uratio;
                % Compute the flame shape and area
                [CI.FM.NL.Model4.q_ratio,CI.FM.NL.Model4.xi,CI.FM.NL.Model4.bashforth_data ] = Fcn_TD_Gequ_interface...
                ( CI.FM.NL.Model4.SU, CI.FM.NL.Model4.xi, CI.FM.NL.Model4.y_vec, CI.TD.dt, 0, CI.FM.NL.Model4.U1, ...
                CI.FM.NL.Model4.area_ratio, CI.TD.uRatio(Var(1):Var(2)),CI.TP.Q * CI.FM.NL.Model4.area_ratio, CI.TP.DeltaHr,CI.FM.NL.Model4.rho1,...
                CI.FM.NL.Model4.bashforth_data,CI.FM.NL.Model4.IT,CI.FM.NL.Model4.time_integration); % In this case CI.TD.uRatio(Var(1):Var(2)) this is a scalar
                CI.FM.NL.Model4.q_ratio_vector(CI.FM.NL.Model4.IT) = CI.FM.NL.Model4.q_ratio;
                CI.TD.qRatio(Var(1):Var(2)) = CI.FM.NL.Model4.q_ratio; % This is a scalar
                
                set(0,'CurrentFigure',CI.FM.NL.Model4.flame_fig)
                plot(CI.FM.NL.Model4.xi,CI.FM.NL.Model4.y_vec)
%                 drawnow
                

            else
                CI.TD.uRatio(Var(1):Var(2)) = Fcn_TD_calculation_uRatio(Var,ss);
                CI.TD.qRatio(Var(1):Var(2)) = Fcn_TD_calculation_qRatio_s(  CI.TD.uRatio,...
                    CI.TD.Green.sys{3},...
                    CI.TD.Green.tauConv2(3),...
                    Var,...
                    CI.TD.taufRem(Var(1):Var(2)),...
                    CI.TD.dt,...
                    CI.TD.Lf);
            end
            % 
            x = CI.TD.IF.Z{ss}*y + CI.TD.IF.Ar*CI.TD.qRatio(Var(1):Var(2)).*CI.TP.RhoCU(ss);
    end
    CI.TD.AP(ss+1,Var(1):Var(2)) = x(1,:);
    CI.TD.AM(ss  ,Var(1):Var(2)) = x(2,:);
end
% -------------------------------------------------------------------------
% --------------------- Right boundary ------------------------------------
% -------------------------------------------------------------------------
% AN-(t) = conv(R2,AN+(t-tauN+))
%
CI.TD.AM(end, Var(1):Var(2))...
        = Fcn_lsim_varied_td(   CI.TD.AP(end,:),...
                                CI.TD.Green.sys{2},...
                                CI.TD.Green.tauConv2(2),...
                                Var,...
                                CI.TP.tau_plus(end) + CI.BC.tau_d2,...
                                CI.TD.dt);
%
% -----------------------------end-----------------------------------------
