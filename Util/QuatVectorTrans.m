function x_global = QuatVectorTrans(quat, x_local)
% Quaternion vector transformation
q_x_local = quaternion(single(0), x_local(1), x_local(2), x_local(3)); % A SINGLE HERE
% q_x_local = quaternion(0, x_local(1), x_local(2), x_local(3));
x_global_temp = (quat.conj*q_x_local)*quat;
x_global_temp_arry = compact(x_global_temp);
x_global = x_global_temp_arry(2:4)';