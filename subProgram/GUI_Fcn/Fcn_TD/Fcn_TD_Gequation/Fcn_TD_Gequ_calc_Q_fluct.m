function [ Q_fluct ] = Fcn_TD_Gequ_calc_Q_fluct( y_vect,SU,xi,u_gutter,Ugs)
%Compute the fluctuating heat release from the flame shape

grad_xi = Fcn_TD_Gequ_calc_grad_xi( u_gutter,SU,xi,y_vec );

flame_area = 0;
diff_y_vec = diff(y_vec);
for runner = 1:length(y_vect)-1 % Integral using trapezoidal method
    flame_area = flame_area +  diff_y_vec(runner) * ...
        pi *(y_vect(runner)*sqrt(1.0+grad_xi(runner).^2) + y_vect(runner+1)*sqrt(1.0+grad_xi(runner+1).^2))/2;
end
flame_area = flame_area * 2; % multiply by 2 so that the flame area for top and bottom half of the duct are accounted for
      
Q_fluct = flame_area *Ugs*RH1*DELTAH(FI);

end

