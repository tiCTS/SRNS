function v_next = VelocityPropagate_LLF_RK4(v0, acc0_LU, acc1_LU, lat, h, M, N, g, dt)
    w_e =  single(7.292115);%e-5;  % rad/s
    w_e_LU = [0; w_e*cos(lat)/100000; w_e*sin(lat)/10000];
    
%     v_k1 = g -( 2* skewMatrix(w_e_LU) + skewMatrix( fun_w_el_LU(v0, M, N, lat, h)))*v0 + acc0_LU;
%     v_k2 = g - ( 2*skewMatrix(w_e_LU) + skewMatrix( fun_w_el_LU((v0+0.5*v_k1*dt),M,N,lat,h)))*(v0+0.5*v_k1*dt) + 0.5*(acc0_LU + acc1_LU);
%     v_k3 = g - ( 2*skewMatrix(w_e_LU) + skewMatrix( fun_w_el_LU((v0+0.5*v_k2*dt),M,N,lat,h)))*(v0+0.5*v_k2*dt) + 0.5*(acc0_LU + acc1_LU);
%     v_k4 = g - ( 2*skewMatrix(w_e_LU) + skewMatrix( fun_w_el_LU((v0 + v_k3*dt),M,N,lat,h)))*(v0 + v_k3*dt) + acc1_LU;

    v_k1 = g -( 2* skewMatrix(w_e_LU)*v0 + skewMatrix( fun_w_el_LU(v0, M, N, lat, h))*v0) + acc0_LU;
    v_k2 = g - ( 2*skewMatrix(w_e_LU)*(v0+0.5*v_k1*dt)  + skewMatrix( fun_w_el_LU((v0+0.5*v_k1*dt),M,N,lat,h))*(v0+0.5*v_k1*dt) )+ 0.5*(acc0_LU + acc1_LU);
    v_k3 = g - ( 2*skewMatrix(w_e_LU)*(v0+0.5*v_k2*dt)  + skewMatrix( fun_w_el_LU((v0+0.5*v_k2*dt),M,N,lat,h))*(v0+0.5*v_k2*dt) )+ 0.5*(acc0_LU + acc1_LU);
    v_k4 = g - ( 2*skewMatrix(w_e_LU)*(v0 + v_k3*dt)  + skewMatrix( fun_w_el_LU((v0 + v_k3*dt),M,N,lat,h))*(v0 + v_k3*dt) )+ acc1_LU;

    v_next = v0 + dt/6 * (v_k1 + 2*v_k2 + 2*v_k3 + v_k4);
    