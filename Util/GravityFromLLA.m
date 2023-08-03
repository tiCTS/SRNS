function gamma_LLF = GravityFromLLA(lat, h)
a1 = 9.7803267715;
a2 = 0.0052790414;
a3 = 0.0000232718;
a4 = -0.000003087691089;
a5 = 0.000000004397731;
a6 = 0.000000000000721;
% Lat in rad

g = a1*(1 + a2*sin(lat)^2 + a3*sin(lat)^4) + (a4 + a5*sin(lat)^2)*h + a6*h^2;

gamma_LLF = [0; 0; -g];