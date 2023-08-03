function wE_LU = fun_wE_LU(lat)
    w_e =  7.292115e-5;  % rad/s
    wE_LU = [0; w_e*cos(lat); w_e*sin(lat)];
    