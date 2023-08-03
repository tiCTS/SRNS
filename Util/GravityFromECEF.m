function gamma_e = GravityFromECEF(p_ECEF)
%% Gravity Vector from ECEF coordinate
%parameters
R = 6378137;  % Radius of the Earth, meter
w_e = 7.292115e-5;  % rad/s
GM = 3986005e8;  % m^3s^-2;
J2 = 0.00108263;
J4 =  -2.37091222e-6;
% implement
psi = atan2(p_ECEF(3), sqrt(p_ECEF(1)^2 + p_ECEF(2)^2));
t = sin(psi);
r = norm(p_ECEF);
a1 = -GM/r^2;
c1 = 1+(3/2)*J2*(R/r)^2;%- (15/8)*J4*(R/r)^4;
c2 = -15/2*J2*(R/r)^2;%+ (105/4)*J4*(R/r)^4;
c3 = -(315/8)*J4*(R/r)^4;
d1 = 1+9/2*J2*(R/r)^2;% - (75/8)*J4*(R/r)^4;
d2 = -(15/2)*J2*(R/r)^2 ;%+ (175/4)*J4*(R/r)^4;
d3 = -(315/8)*J4*(R/r)^4;

% gamma_e = a1/r * [(c1+c2*t^2+c3*t^4)*p_ECEF(1);  (c1+c2*t^2+c3*t^4)*p_ECEF(2); (d1+d2*t^2+d3*t^4)*p_ECEF(3)] + ...
%     [w_e^2*p_ECEF(1); w_e^2*p_ECEF(2); 0];   % J2 & J4

gamma_e = a1/r * [(c1+c2*t^2)*p_ECEF(1);  (c1+c2*t^2)*p_ECEF(2); (d1+d2*t^2)*p_ECEF(3)] + ...
    [w_e^2*p_ECEF(1); w_e^2*p_ECEF(2); 0];