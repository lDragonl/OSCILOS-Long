function y = Fcn_calculation_envelope(y0,nPadding,filter)
% This function is used to calculate the envelope of a singal
% y0 is the original signal
% nPadding is the length of zero padding signals before and after the
% original signal y0
% filter is a field, which contains the filter information:
% fs fPass, fStop, Rp, Rs, filterType
% for detail information, please refer to the function 
%
% last edit: 2014.11.14 14:21
%
IsRowArray = 1;
SizeY = size(y0);
if SizeY(1)>SizeY(2)
    y0 = y0';               % The array should be 1xn style
    IsRowArray = 0;
end
%
% -------------------------------------------------------------------------
% The third input should be the filter information, 
% filter.fs           = 1e5;
% filter.fPass        = 300;
% filter.fStop        = 600;
% filter.Rp           = 1;
% filter.Rs           = 40;
% filter.filterType   = 1;
%
yP              = [zeros(1,nPadding), y0, zeros(1,nPadding)]; 

                                           
yEnv    = abs(hilbert(detrend(yP)));
% -------------------------------------------------------------------------
y       = yEnv(nPadding+1:end-nPadding);
%
% ----------------------------end------------------------------------------



