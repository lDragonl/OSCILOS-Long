function [ grad_xi ] = Fcn_TD_Gequ_calc_grad_xi( u_gutter,SU,xi,y_vec )
%Simple backward euler first oder graidnt method. The boundary condition is
%applied such that the flame remains anchored at y(1) if SU < u_gutter. If
%this is not the case, then the gradinet at y(1) is set to zero

% Compute BC
if u_gutter >= SU
   BC_grad = sqrt(( u_gutter/SU)^2 - 1);
end
    
grad_xi = [BC_grad, diff(xi)./diff(y_vec)];

end

