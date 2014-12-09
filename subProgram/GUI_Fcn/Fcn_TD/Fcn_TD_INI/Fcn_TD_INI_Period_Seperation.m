function Fcn_TD_INI_Period_Seperation(nPeriod,RatioGapPadding)
% In case the flame describing function is determined by the envelope of
% velocity ratio, uRatio, iteration-convergence is necessary. 
% The whole calculation can be divided to several intervals with
% equal length. 
% a gap is necessary between two neighbouring intervals
% last edit: 2014-11-13
%
if nargin == 0
    % default value
    nPeriod = 10;
    RatioGapPadding = 0;
end
global CI
if CI.FM.NL.style == 4
    CI.TD.nPeriod               = CI.TD.nRound;  
else
    CI.TD.nPeriod               = nPeriod;                                          % number of interval
end

    CI.TD.numGap                = round(CI.TD.nRound./CI.TD.nPeriod);               % number of Gap counting
switch CI.FM.NL.style                                                           % In case the envelope of uRatio is necessary, a padding is used. 
    case 3
        CI.TD.RatioGapPadding  = RatioGapPadding;                               % ratio of padding length  
    otherwise
        CI.TD.RatioGapPadding  = 0;
end
CI.TD.numGapPadding         = round(CI.TD.RatioGapPadding.*CI.TD.numGap);       % number of padding
%
for ss = 1:CI.TD.nPeriod
    CI.TD.indexPeriod(ss,1) = (ss-1).*CI.TD.numGap + 1;                         % start of one period
    CI.TD.indexPeriod(ss,2) = ss.*CI.TD.numGap + CI.TD.numGapPadding;           % end of one period
    if CI.TD.indexPeriod(ss,2) > CI.TD.nRound
        CI.TD.indexPeriod(ss,2) = CI.TD.nRound;
    end
end                                 
assignin('base','CI',CI);
% ----------------------------end------------------------------------------
