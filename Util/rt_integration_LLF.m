%% Rotation matrix integration in LLF
% This script integrate the rotation matrix in a period of time (dt) with
% given angular velocity (w) and initial rotation matrix (R) in LLF
% The rotation matrix R transform the vector from body frame to LLF.
% Input:
%   R : initial rotation matrix ,3x3 matrix
%   w : angular velocity measurement from gyro, 3x1 vector
%   wE_l : Earth angular velocity in LLF, 3x1 vector
%   w_el_l : LLF frame angular velocity w.r.t. E frame, 3x1 vector
%   dt : time period
% Output:
%   R_next : the rotation matirx after dt, 3x3 matrix
function R_next = rt_integration_LLF(R, w, wE_l, w_el_l,  dt)
    S_ib_b = skewMatrix(dt.*w);
    w_il_b = R*(wE_l + w_el_l);
    S_il_b = skewMatrix( dt.* w_il_b);

    % the coefficients
    s_ib = rt_exp_sin(dt.*w);
    c_ib = rt_exp_cos(dt.*w);

    s_il = rt_exp_sin(dt.*w_il_b);
    c_il = rt_exp_cos(dt.*w_il_b);

    R_next = R*(eye(3) + s_ib*S_ib_b + c_ib*S_il_b^2)/(eye(3) + s_il*S_il_b + c_il*S_il_b^2);
end