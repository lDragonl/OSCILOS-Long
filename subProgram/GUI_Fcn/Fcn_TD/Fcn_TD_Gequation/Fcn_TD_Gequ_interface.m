function [ q_ratio,xi ] = Fcn_TD_Gequ_interface( SU, xi, y_vec, dt, current_time, U1, area_ratio, uratio,Q_mean,delta_h,rho1 )
%This is the G-equation interface with OSCILOS. It asks in all of the
%required variables, and gives out the heat release ratio to be used in
%the acoustics calculations

% The G_equation is a three step operation, First findf the velocity right before the flame
u_gutter = Fcn_TD_Gequ_calc_ugutter( U1,area_ratio,uratio,current_time);

% Then compute the new flame shape
% For this we have a choice of different numerical methods
time_integration = 1;
switch time_integration 
    case(1)
        [ dxidt ] = Fcn_TD_Gequ_calc_dxidt( u_gutter,SU,xi,y_vec );
        [ xi ] = Fcn_TD_Gequ_forward_euler( xi,dt,dxidt );
    case(2)
        [ xi ] = Fcn_TD_Gequ_RK(SU,xi,y_vec,dt,current_time,U1,area_ratio,uratio,u_gutter );
end

% Then compute the heat release from the area
[ Q_fluct ] = Fcn_TD_Gequ_calc_Q_fluct( y_vec,SU,xi,u_gutter, Q_mean, delta_h,rho1,y_vec);

% and finally compute the Q_ratio
q_ratio = Q_fluct/Q_mean;

end

