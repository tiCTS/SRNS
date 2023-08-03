function v_next = VelocityPropagate_LLF_ZOH(v0, acc0_LU, lat, h, M, N, g, dt)
    w_e =  single(7.292115);%e-5;  % rad/s
    w_e_LU = [0; w_e*cos(lat)/100000; w_e*sin(lat)/10000];

    dv = g -( 2* skewMatrix(w_e_LU)*v0 + skewMatrix( fun_w_el_LU(v0, M, N, lat, h))*v0) + acc0_LU;

    v_next = v0 + dt * dv;