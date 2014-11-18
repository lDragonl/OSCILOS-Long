function [ xi_steady ] = Fcn_TD_Gequ_steady_flame( U_gutter_steady,SU,y_vec )
%This function computes the one sided steady flame shape of the flame using the vale
%of the laminar burning velocity SU and the mean flow. U1

slope = sqrt((U_gutter_steady/SU)^2 - 1);
xi_steady = y_vec * slope - y_vec(1) * slope;

end

