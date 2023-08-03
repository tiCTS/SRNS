function M = M_LLF(lat)
    a_Earth = 6378137; % WGS84 Earth semi-major axis (m).
    e_sq_Earth = 6.694379990e-3; % WGS84 Earth eccentricity squared.

    M = (a_Earth*(1-e_sq_Earth))/(1 - e_sq_Earth*sin(lat)^2)^(3/2);