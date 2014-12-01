function [ grad_xi ] = Fcn_TD_Gequ_calc_grad_xi( u_gutter,SU,xi,y_vec )
%Simple backward euler first oder graidnt method. The boundary condition is
%applied such that the flame remains anchored at y(1) if SU < u_gutter. If
%this is not the case, then the gradinet at y(1) is set to zero

% Compute BC
if u_gutter >= SU && xi(:,1)
   BC_grad = (xi(:,1) >= 0 ) .* sqrt(( u_gutter/SU).^2 - 1); % Implicit test is to determine whether or not we have reached the anchoring point
else
   BC_grad = zeros(size(u_gutter));
end

prod_vec = ones(length(u_gutter),1);
grad_xi = cat(2,BC_grad, diff(xi,1,2)./(prod_vec * diff(y_vec)));

end

