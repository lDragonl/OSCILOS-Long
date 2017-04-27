
U_mean = 9.15;
Lf = 0.1;

beta = 5.660158275081317;
gamma = 0.5;
%SL = U_mean/sqrt(1+beta^2/gamma^2);
K = 0.8;

s = linspace(0.1,1000,1001)*2*pi*1i;
eta = K*beta^2/(beta^2+gamma^2);
St = imag(s)*Lf/U_mean;
%St = s*CI.FM.HP{indexHP}.GEQU_CONV.Lf/(1i*CI.FM.HP{indexHP}.GEQU_CONV.Ugs);
St2 = St*(beta^2+gamma^2)/beta^2;

F = 2*(exp(1i*eta*St2).*(gamma-1i*eta*St2) - eta*exp(1i*St2).*(gamma-1i*St2) + gamma*(eta-1))./...
    (eta*St2.^2*(2-gamma)*(eta-1));

plot(s/(2*pi*1i),abs(F))