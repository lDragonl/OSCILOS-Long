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
    if ~isempty(CI.CD.indexHP)    % if there are heat perturbations
        
        if CI.FM.indexFM(1)== 4 % G-equation case, at the moment all flame models must be the same for all flames
            % AO: This is because this G-equation code has been writte in
            % terms of "flame vectors". It would be better to code it for
            % single flames and then loop through them.
            CI.TD.IT = nn;
            %
            % AO: G-equation convective case
        elseif CI.FM.indexFM(1) == 5
            CI.TD.IT = nn;            
        end
        
    end
    
    Fcn_TD_calculation_one_gap(Var)
    waitbar(nn/CI.TD.nRound);
    drawnow
    
end

close(hWaitBar)
assignin('base','CI',CI);

%
% ----------------------------end------------------------------------------