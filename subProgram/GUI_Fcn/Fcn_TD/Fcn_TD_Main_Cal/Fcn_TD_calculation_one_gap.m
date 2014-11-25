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
CI.TD.AP(1, Var(1):Var(2)) = CI.TD.AP(1, Var(1):Var(2)) + CI.TD.pNoiseBG(Var(1):Var(2));
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
            CI.TD.uRatio(Var(1):Var(2)) = Fcn_TD_calculation_uRatio(Var,ss);
            if CI.FM.NL.style == 4
            	CI.TD.qRatio(Var(1):Var(2)) = Fcn_TD_calculation_qRatio_s(  CI.TD.uRatio,...
            else
                CI.TD.qRatio(Var(1):Var(2)) = Fcn_TD_calculation_qRatio(CI.TD.uRatio,CI.TD.Green.y{3},Var,CI.TD.dt,CI.TD.Lf);
            end                                                            CI.TD.Green.sys{3},...
                                                                        CI.TD.Green.tauConv2(3),...
                                                                        Var,...
                                                                        CI.TD.taufRem(Var(1):Var(2)),...
                                                                        CI.TD.dt,...
                                                                        CI.TD.Lf);
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
