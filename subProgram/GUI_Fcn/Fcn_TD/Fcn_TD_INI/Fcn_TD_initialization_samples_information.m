function Fcn_TD_initialization_samples_information(fs,tEndRaw,uRatioMax,nGapRatio)
% This function is used to initialize the samples information
% The sampling rate, time end, nGap, and some time samples are initialize
% in this sub-code. Its parent program is "Fcn_TD_initialization".
% fs, tEndRaw,nGapRatio,uRatioMax
% are the input. More input can be added in the future version.
%
% last edited: 2014-11-12 16:41 
%
global CI
% -------------------------------------------------------------------------
%
CI.TD.fs                = fs;                                               % Sampling frequency
CI.TD.tEndRaw           = tEndRaw;                                          % a raw time end, which may be slightly changed
                                                                            % in the following code
%
% -------------------------------------------------------------------------
%
CI.TD.dt                = 1/CI.TD.fs;                                       % time step
CI.TD.NRaw              = round(CI.TD.tEndRaw*CI.TD.fs);                    % raw number of samples after time "0".
%
% -------------------------------------------------------------------------

[CI.TD.tauMin, CI.TD.taufMin] = Fcn_TD_minimum_time_delay(uRatioMax);       % CI.TD.taufMin is the minimum time delay of the FTF, 
                                                                            % which becomes very important when the time delay 
                                                                            % of the flame describing functing function varies
                                                                            % with velocity ratio
                                                                            % CI.TD.taufMin is a minimum time delay to judge 
                                                                            % the maximum number of samples in one calculation period.

%
CI.TD.nGapMax           = floor(CI.TD.tauMin.*CI.TD.fs);                    % nGapMax is the maximum time steps which can be calculated in one loop
CI.TD.nGap              = round(nGapRatio.*CI.TD.nGapMax);                  % This value should be smaller than the maximum value
CI.TD.nRound            = ceil(CI.TD.NRaw./CI.TD.nGap);                     % This is the number of loops to be calculated
CI.TD.nCal              = CI.TD.nRound*CI.TD.nGap;                          % This is the final sampling number after time zero
CI.TD.tSp               = (1:CI.TD.nRound*CI.TD.nGap).*CI.TD.dt;            % This is final time samples afte time zero
% ------------------
CI.TD.tauPadding        = Fcn_TD_Padding_time_delay;                        % accounting for the time delay, zero padding before the time zero is necessary
CI.TD.nPadding          = round(CI.TD.tauPadding.*CI.TD.fs);                % number of padding signal
CI.TD.tSpPadding        = (-CI.TD.nPadding+1:0).*CI.TD.dt;                  % padding time samples
% ------------------
CI.TD.tSpTotal          = cat(2,CI.TD.tSpPadding,CI.TD.tSp);                % Total time samples
CI.TD.nTotal            = CI.TD.nCal + CI.TD.nPadding;                      % number of total time samples
%
assignin('base','CI',CI);
%
% -----------------------------end-----------------------------------------