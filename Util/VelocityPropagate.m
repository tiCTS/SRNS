function v_next = VelocityPropagate(v0, p_ECEF_0, f0, f1, dt)
% Velocity integration using RK4
w_e = 7.292115e-5;  % rad/s
w_IE_E = [0;0;w_e];  % rad/s
vk1 = -skewMatrix(w_IE_E)*skewMatrix(w_IE_E)*p_ECEF_0 - 2*skewMatrix(w_IE_E)*v0 + f0;
vk2 = -skewMatrix(w_IE_E)*skewMatrix(w_IE_E)*p_ECEF_0 - 2*skewMatrix(w_IE_E)*(v0 + 0.5*dt*vk1) + 0.5*(f0 + f1);
vk3 = -skewMatrix(w_IE_E)*skewMatrix(w_IE_E)*p_ECEF_0 - 2*skewMatrix(w_IE_E)*(v0 + 0.5*dt*vk2) + 0.5*(f0 + f1);
vk4 = -skewMatrix(w_IE_E)*skewMatrix(w_IE_E)*p_ECEF_0 - 2*skewMatrix(w_IE_E)*(v0 + dt*vk3) + f1;
v_next = v0 + dt/6 * (vk1 + 2*vk2 + 2*vk3 + vk4);