function Fcn_TD_INI_background_noise(NoiseInfo)
% -------------------------------------------------------------------------
% Background noise
% This function is used set the background noise, which can be used to 
% trigger the time domain simulations 
% The pressure is defaultly added at the inlet boundary 
%
% first created: 2014-11-12 17:14
% last edited: 2014-11-21
%
% check input 
global CI
if nargin == 0
    % default value
    NoiseInfo.level     = -40;                                                  % pNoise level 0.01 Pa
    NoiseInfo.t(1)      = CI.TD.tSpTotal(1);
    NoiseInfo.t(end)    = CI.TD.tSpTotal(end);
else
    if NoiseInfo.t(1) < CI.TD.tSpTotal(1)
        NoiseInfo.t(1) = CI.TD.tSpTotal(1);
    end
    if NoiseInfo.t(end) > CI.TD.tSpTotal(end)
        NoiseInfo.t(end) = CI.TD.tSpTotal(end);
    end
end
CI.TD.NoiseInfo         = NoiseInfo;                                            % background noise level = 10^(value/20) Pa;
CI.TD.pNoiseBG          = wgn(1,CI.TD.nTotal,CI.TD.NoiseInfo.level);            % background noise
N1                      = round((NoiseInfo.t(1)-CI.TD.tSpTotal(1)+1).*CI.TD.fs);
N2                      = round((NoiseInfo.t(end)-CI.TD.tSpTotal(1)+1).*CI.TD.fs);
CI.TD.pNoiseBG(1:N1)    = 0;
CI.TD.pNoiseBG(N2:end)  = 0;
%
assignin('base','CI',CI);
%
% ------------------------end----------------------------------------------