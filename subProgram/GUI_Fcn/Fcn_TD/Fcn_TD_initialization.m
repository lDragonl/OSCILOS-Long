function Fcn_TD_initialization(fs,tEndRaw,uRatioMax,nGapRatio,nPeriod,RatioGapPadding,NoiseLevel)
global CI
% This function is used to initilize the time domain calculation code
% The sampling frequency and the raw time ends,
% and ...
% last edit: 2014-11-13
%
% -------------------------------------------------------------------------
% check input 
if nargin == 0
    % default value
    fs              = 1e5;
    tEndRaw         = 1;
    uRatioMax       = 1;
    nGapRatio       = 1;
    nPeriod         = 10;
    RatioGapPadding = 9;
    NoiseLevel      = -40;
end
%
addpath(genpath('./'))                                                      % add directories to search path
% -------------------------------------------------------------------------
%
CI.TD.uRatioMax = uRatioMax;                                                % Set an uplimit of uRatio, which ensures that
%                                                                           % the time domain simulation results are not unreasonable. 
%
Fcn_TD_initialization_matrix_preprocessing                                  % preprocessing some matrix connecting the components
%                                                                           % at two sides of an interface
%
Fcn_TD_initialization_samples_information(fs,tEndRaw,uRatioMax,nGapRatio)   % initialize some samples
%       
Fcn_TD_initialization_waves                                                 % initialize waves
%
Fcn_TD_initialization_flame_model                                           % some flame describing function model information
%
Fcn_TD_initialization_background_noise(NoiseLevel)                          % background noise information
%
Fcn_TD_initialization_Green_function                                        % Green's function information
% 
Fcn_TD_initialization_constant_variables                                    % pre-calculate some constant values to improve the calculation speed
%
Fcn_TD_initialization_Period_Seperation(nPeriod,RatioGapPadding)            %
%
assignin('base','CI',CI);
%
% ------------------------ end --------------------------------------------