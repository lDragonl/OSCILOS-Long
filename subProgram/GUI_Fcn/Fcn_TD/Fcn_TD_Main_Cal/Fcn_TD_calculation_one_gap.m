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
if CI.TD.ExtForceInfo.indexPos == 1
    CI.TD.AP(1, Var(1):Var(2)) = CI.TD.AP(1, Var(1):Var(2)) + CI.TD.pExtForce(Var(1):Var(2));  % external driving;
end
%
% -------------------------------------------------------------------------
% ------------------ interfaces between two sections ----------------------
% -------------------------------------------------------------------------
indexHA_num = 0; % Initilise the index which counts the number of sections which have either mean heat release or fluctuating heat release
indexHP_num = 0; % Initilise the index which counts the number of sections which have fluctuating heat release only
indexHPTF_num = 0; % Initilise the index which counts the number of sections of which the flame model is expressed as flame transfer function

for ss = 1:CI.TP.numSection-1
    y(1,1:CI.TD.nGap)   = Fcn_interp1_varied_td(CI.TD.AP(ss,:),    Var, CI.TP.tau_plus(ss),    CI.TD.dt);
    y(2,1:CI.TD.nGap)   = Fcn_interp1_varied_td(CI.TD.AM(ss+1,:),  Var, CI.TP.tau_minus(ss+1), CI.TD.dt);
    y(3,1:CI.TD.nGap)   = 0;
    if CI.TD.ExtForceInfo.indexPos == 2*ss % the loudspeaker is located at the left side of interface
        y(1,1:CI.TD.nGap) = y(1,1:CI.TD.nGap) + 0.5*CI.TD.pExtForce(Var(1):Var(2));
    end
    if CI.TD.ExtForceInfo.indexPos == 2*ss+1 % the loudspeaker is located at the left side of interface
        y(2,1:CI.TD.nGap) = y(2,1:CI.TD.nGap) + 0.5*CI.TD.pExtForce(Var(1):Var(2));
    end
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
            indexHP_num = indexHP_num + 1;
            
            %AO: Added for the convective case
            if CI.FM.indexFM(indexHP_num) == 5  % In the case of the G-equation, the value of qratio is computed from a pde, not from a transfer function using Greens function
                
                
                CI.TD.uRatio(indexHP_num,Var(1):Var(2)) = Fcn_TD_calculation_uRatio(Var,ss, indexHP_num);
                %disp(num2str(u_ac));
                
                isReINI_Global = 0;
                if rem(CI.TD.IT, 10) == 0
                    isReINI_Global = 1;
                end
                
                [CI.FM.HP{indexHP_num}.GEQU_CONV.dataGhost, CI.FM.HP{indexHP_num}.GEQU_CONV.FF, CI.FM.HP{indexHP_num}.GEQU_CONV.vyPrime_no_stencil] = ...
                    FcnLLS_MainCal_Dynamic(CI.TD.dt*(CI.TD.IT-1), CI.FM.HP{indexHP_num}.GEQU_CONV.SET, CI.FM.HP{indexHP_num}.GEQU_CONV.dataGhost, CI.FM.HP{indexHP_num}.GEQU_CONV.gridGhost, CI.FM.HP{indexHP_num}.GEQU_CONV.vyPrime_no_stencil, CI.TD.uRatio(indexHP_num,Var(1):Var(2)), isReINI_Global);
                
                
%                 [ CI.FM.HP{indexHP_num}.GEQU_CONV.q_ratio, CI.FM.HP{indexHP_num}.GEQU_CONV.xi, CI.FM.HP{indexHP_num}.GEQU_CONV.uflow ] = ...
%                     Fcn_TD_Conv_Gequ_interface( CI.FM.HP{indexHP_num}.GEQU_CONV.SU, CI.FM.HP{indexHP_num}.GEQU_CONV.xi, ...
%                     CI.FM.HP{indexHP_num}.GEQU_CONV.y_vec, CI.TD.dt, CI.TD.dt*(CI.TD.IT-1), CI.FM.HP{indexHP_num}.GEQU_CONV.U1, ...
%                     CI.FM.HP{indexHP_num}.GEQU_CONV.area_ratio, CI.FM.HP{indexHP_num}.GEQU_CONV.uflow, CI.TP.Q(indexHA_num), ...
%                     CI.TP.DeltaHr(indexHA_num), CI.FM.HP{indexHP_num}.GEQU_CONV.rho1 , CI.FM.HP{indexHP_num}.GEQU_CONV.x_flow, ...
%                     CI.FM.HP{indexHP_num}.GEQU_CONV.UC, 0, CI.TD.uRatio(indexHP_num,Var(1):Var(2)), 0);
                %CI.FM.HP{indexHP_num}.GEQU_CONV.UC, 1, 0.2, 150);

                CI.FM.HP{indexHP_num}.GEQU_CONV.saved_rs(CI.TD.IT) = CI.FM.HP{indexHP_num}.GEQU_CONV.FF.xFFInterp;
                CI.FM.HP{indexHP_num}.GEQU_CONV.saved_xs(CI.TD.IT) = CI.FM.HP{indexHP_num}.GEQU_CONV.FF.yFFInterp;
                if(CI.TD.IT ==1)
                    CI.TD.qZer0 = CI.FM.HP{indexHP_num}.GEQU_CONV.FF.A;
                end
                CI.TD.qRatio(indexHP_num ,Var(1):Var(2)) = (CI.FM.HP{indexHP_num}.GEQU_CONV.FF.A-CI.TD.qZer0)/CI.TD.qZer0;
                
                
            elseif CI.FM.indexFM(indexHP_num) == 4  % In the case of the G-equation, the value of qratio is computed from a pde, not from a transfer function using Greens function
                
                CI.TD.uRatio(indexHP_num,Var(1):Var(2)) = Fcn_TD_calculation_uRatio(Var,ss, indexHP_num); % compute the values of uratio in this section
                
                % Compute the flame shape and area
                [CI.FM.HP{indexHP_num}.GEQU.q_ratio,CI.FM.HP{indexHP_num}.GEQU.xi,CI.FM.HP{indexHP_num}.GEQU.bashforth_data ] = Fcn_TD_Gequ_interface...
                    ( CI.FM.HP{indexHP_num}.GEQU.SU, CI.FM.HP{indexHP_num}.GEQU.xi, CI.FM.HP{indexHP_num}.GEQU.y_vec, CI.TD.dt, 0, CI.FM.HP{indexHP_num}.GEQU.U1, ...
                    CI.FM.HP{indexHP_num}.GEQU.area_ratio, CI.TD.uRatio(indexHP_num,Var(1):Var(2)),CI.TP.Q(indexHA_num), CI.TP.DeltaHr(indexHA_num),...
                    CI.FM.HP{indexHP_num}.GEQU.rho1,CI.FM.HP{indexHP_num}.GEQU.bashforth_data,CI.TD.IT,CI.FM.HP{indexHP_num}.GEQU.time_integration); % In this case CI.TD.uRatio(indexHA_num,Var(1):Var(2)) is a scalar
                CI.FM.HP{indexHP_num}.GEQU.saved_xi(CI.TD.IT,:) = CI.FM.HP{indexHP_num}.GEQU.xi;
                
                
                
                CI.TD.qRatio(indexHP_num ,Var(1):Var(2)) = CI.FM.HP{indexHP_num}.GEQU.q_ratio; % This is a scalar
                
            elseif CI.FM.indexFM(indexHP_num) < 3
                indexHPTF_num = indexHPTF_num + 1;
                CI.TD.uRatio(indexHP_num,Var(1):Var(2)) = Fcn_TD_calculation_uRatio(Var,ss,indexHP_num);
                CI.TD.qRatio(indexHP_num,Var(1):Var(2)) = Fcn_TD_calculation_qRatio_s(  CI.TD.uRatio(indexHP_num,:),...
                    CI.TD.Green.sys{indexHPTF_num + 2},...
                    CI.TD.Green.tauConv2(indexHPTF_num + 2),...
                    Var,...
                    CI.TD.taufRem(indexHP_num ,Var(1):Var(2)),...
                    CI.TD.dt,...
                    CI.TD.Lf(indexHP_num,:),...
                    indexHP_num);
            end
            %
            
            x = CI.TD.IF.Z{ss}*y + CI.TD.IF.Ar{indexHP_num}*CI.TD.qRatio(Var(1):Var(2)).*CI.TP.RhoCU(ss);
            
            if CI.TD.ExtForceInfo.indexPos == 2*ss % the loudspeaker is located at the left side of interface
                x(2,1:CI.TD.nGap) = x(2,1:CI.TD.nGap) + 0.5*CI.TD.pExtForce(Var(1):Var(2));
            end
            if CI.TD.ExtForceInfo.indexPos == 2*ss+1 % the loudspeaker is located at the left side of interface
                x(1,1:CI.TD.nGap) = x(1,1:CI.TD.nGap) + 0.5*CI.TD.pExtForce(Var(1):Var(2));
            end
            %AO: checks that every section is assigned to a certain value
        otherwise
            
            error('Section type undefined')
            
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
CI.TD.AM(end, Var(1):Var(2)) = CI.TD.AM(end, Var(1):Var(2));% + CI.TD.pNoiseBG(Var(1):Var(2)); % add additional noise to wave propagating in direction of flow
if CI.TD.ExtForceInfo.indexPos == 2*(length(CI.CD.x_sample) - 1)
    CI.TD.AM(end, Var(1):Var(2)) = CI.TD.AM(end, Var(1):Var(2))+ CI.TD.pExtForce(Var(1):Var(2));  % external driving;
end

% AO: debugging purposes
if(mod(CI.TD.IT,100)==0)
    figure(1)
    clf
    subplot(3,1,1)
    plot(CI.TD.tSpTotal,CI.TD.qRatio,'k-')
    subplot(3,1,2)
    hold on
    xs = CI.FM.HP{1}.GEQU_CONV.saved_xs(CI.TD.IT);
    rs = CI.FM.HP{1}.GEQU_CONV.saved_rs(CI.TD.IT);
    for ii=1:length(xs)
        plot(xs{ii},rs{ii},'r-')
    end
    subplot(3,1,3)
    plot(CI.FM.HP{1}.GEQU_CONV.vyPrime_no_stencil,'k-')
end
assignin('base','CI',CI);
%
% -----------------------------end-----------------------------------------
