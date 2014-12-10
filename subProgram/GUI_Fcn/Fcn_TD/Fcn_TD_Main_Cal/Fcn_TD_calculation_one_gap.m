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
% AP is the pressure wave propagating in the direction of flow
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
indexHA_num = 0; % Initilise the index which counts the number of sections which have either mean heat release or fluctuating heat release
indexHF_num = 0; % Initilise the index which counts the number of sections which have fluctuating heat release only
for ss = 1:CI.TP.numSection-1 
    y(1,1:CI.TD.nGap)   = Fcn_interp1_varied_td(CI.TD.AP(ss,:),    Var, CI.TP.tau_plus(ss),    CI.TD.dt);
    y(2,1:CI.TD.nGap)   = Fcn_interp1_varied_td(CI.TD.AM(ss+1,:),  Var, CI.TP.tau_minus(ss+1), CI.TD.dt);
    y(3,1:CI.TD.nGap)  = 0; 
    switch CI.CD.SectionIndex(ss+1)
        case 0   % only area change
            % ---------------------------
            %     [ A2+ ]       [ A1+ ]
            %     [ A1- ] =  [Z][ A2- ]
            %     [ E2  ]       [ E1  ]
            %
            x = CI.TD.IF.Z{ss}*y;
            
        case 10   % Mean heat release. 
                  % The fluctuating pressure waves remain unaffected, but
                  % the counter indexHA_num is increased
            % ---------------------------
            %     [ A2+ ]       [ A1+ ]
            %     [ A1- ] =  [Z][ A2- ]
            %     [ E2  ]       [ E1  ]
            %
            indexHA_num = indexHA_num + 1;
            x = CI.TD.IF.Z{ss}*y;    
            
            % ---------------------------
        case 11   % with heat addition and perturbations, the counter indexHA_num is increased
            % ---------------------------
            %     [ A2+ ]       [ A1+ ]   [    ]
            %     [ A1- ] =  [Z][ A2- ] + [ Ar ] conv(FTF,u*rho*c)
            %     [ E2  ]       [ E1  ]   [    ]
            %     u*(rho*c) = A1+(t-tauf) - A1-(t-tauf)
            %
            indexHA_num = indexHA_num + 1;
            indexHF_num = indexHF_num + 1;
            if CI.FM.NL.style == 4 % In the case of the G-equation, the value of qratio is computed from a pde, not from a transfer function using Greens function
                
                CI.TD.uRatio(Var(1):Var(2)) = Fcn_TD_calculation_uRatio(Var,ss); % compute the values of uratio in every section
                
                % Compute the flame shape and area
                [CI.FM.NL.Model4.q_ratio(indexHF_num),CI.FM.NL.Model4.xi(indexHF_num,:),CI.FM.NL.Model4.bashforth_data(indexHF_num,:,:) ] = Fcn_TD_Gequ_interface...
                ( CI.FM.NL.Model4.SU(indexHF_num), CI.FM.NL.Model4.xi(indexHF_num,:), CI.FM.NL.Model4.y_vec(indexHF_num,:), CI.TD.dt, 0, CI.FM.NL.Model4.U1(indexHF_num), ...
                CI.FM.NL.Model4.area_ratio(indexHF_num), CI.TD.uRatio(indexHA_num,Var(1):Var(2)),CI.TP.Q(indexHA_num) * CI.FM.NL.Model4.area_ratio(indexHF_num), CI.TP.DeltaHr(indexHA_num),...
                CI.FM.NL.Model4.rho1(indexHF_num),CI.FM.NL.Model4.bashforth_data(indexHF_num,:,:),CI.FM.NL.Model4.IT,CI.FM.NL.Model4.time_integration); % In this case CI.TD.uRatio(indexHA_num,Var(1):Var(2)) this is a scalar
                
                CI.TD.qRatio(indexHF_num,Var(1):Var(2)) = CI.FM.NL.Model4.q_ratio(indexHF_num); % This is a scalar
              
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
