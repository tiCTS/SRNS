function q_next = Quaternion1stInt (w_t1, w_t2, quat, dt)
%% Quaternion 1st-order integration
w_dot = (w_t2 - w_t1)/dt;
w_bar = w_t1 + 0.5*w_dot*dt;
% incremental rotation with w_bar
q_inc = [ cos(0.5*norm(w_bar)*dt)  (1/norm(w_bar))*sin(0.5*norm(w_bar)*dt)*w_bar'];

q_next = quat * quaternion( q_inc + dt^2/24 * [0; cross(w_t1, w_t2)]');