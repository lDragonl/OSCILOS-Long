function Fcn_TD_INI_matrix_preprocessing
% This function is used to calculate the matrix connecting the components
% at two sides of an interface. The inputs are the known values and the
% outputs are the values to be calculated. The matrix in the technical
% report should be reorganized.
% 
% last edit: 2014-11-13
%
global CI
for ss = 1:CI.TP.numSection-1 
    switch CI.CD.index(ss+1)
        case {0,10} % case of no heat release, or mean heat release only
            CI.TD.IF.Z{ss} = Fcn_TD_matrix_reorganization_no_flame(CI.TPM.BC{ss}); 
        case {11} % case of heat release perturbations
            % heat addition coefficient
            K       = CI.TP.DeltaHr./CI.TP.c_mean(2,ss+1)./CI.TP.c_mean(1,ss)./CI.TP.Theta(ss);
            % right side, first term
            B1a     = CI.TPM.B1{2,ss}*(CI.TPM.B2{1,ss}\CI.TPM.B1{1,ss});
            BC1a    = B1a*CI.TPM.C1;
            % left side
            BC2     = CI.TPM.B2{2,ss}*CI.TPM.C2;
            [CI.TD.IF.Z{ss},CI.TD.IF.Ar]...
                = Fcn_TD_matrix_reorganization_with_flame(BC1a,BC2,K);
    end
end
assignin('base','CI',CI);
% -------------------------------------------------------------------------
%
function M2 = Fcn_TD_matrix_reorganization_no_flame(M1)
% change the matrix M1
%     [ A2+ ]       [ A1+ ]
%     [ A2- ] = [M1][ A1- ]
%     [ E2  ]       [ E1  ]
% to:
%     [ A2+ ]       [ A1+ ]
%     [ A1- ] = [M2][ A2- ]
%     [ E2  ]       [ E1  ]
Z       = M1;
ZLeft   = [     1   -Z(1,2)     0;...
                0   -Z(2,2)     0;...
                0   -Z(3,2)     1];
ZRight  = [     Z(1,1)      0       Z(1,3);...
                Z(2,1)      -1      Z(2,3);...
                Z(3,1)      0       Z(3,3)];
M2      = ZLeft\ZRight;
% 
% -------------------------------------------------------------------------
%
function [M3,Ar] = Fcn_TD_matrix_reorganization_with_flame(M1,M2,K)
% change the matrix 
%           
%         [ A2+ ]       [ A1+ ]   [ 0 ]
%     [M2][ A2- ] = [M1][ A1- ] + [ 0 ] FTF*u
%         [ E2  ]       [ E1  ]   [ K ]
% to: 
%         [ A2+ ]       [ A1+ ]   [    ]
%         [ A1- ] = [M3][ A2- ] + [ Ar ] FTF*u
%         [ E2  ]       [ E1  ]   [    ]
% where Ar is an array
Ar1     = [0 0 K]';
Z       = M2\M1;
Ar2     = M2\Ar1;
ZLeft   = [     1   -Z(1,2)     0;...
                0   -Z(2,2)     0;...
                0   -Z(3,2)     1];
ZRight  = [     Z(1,1)      0       Z(1,3);...
                Z(2,1)      -1      Z(2,3);...
                Z(3,1)      0       Z(3,3)];
M3 = ZLeft\ZRight;
Ar = ZLeft\Ar2;
% 
% ----------------------------end-------------------------------------------