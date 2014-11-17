function Fcn_TD_initialization_background_noise(NoiseLevel)
% -------------------------------------------------------------------------
% Background noise
% This function is used set the background noise, which can be used to 
% trigger the time domain simulations 
% The pressure is defaultly added at the inlet boundary 
%
% last edited: 2014-11-12 17:14
%
% check input 
if nargin == 0
    % default value
    NoiseLevel = -40;                                                       % pNoise level 0.01 Pa
end
global CI
CI.TD.NoiseLevel          = NoiseLevel;                                     % background noise level = 10^(value/20) Pa;
CI.TD.pNoiseBG            = wgn(1,CI.TD.nTotal,CI.TD.NoiseLevel);           % background noise
%
%
assignin('base','CI',CI);
%
% ------------------------end----------------------------------------------