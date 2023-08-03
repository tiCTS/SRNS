%% RT_Expansion_Cos
% This function computes the coefficient of cos term of Taylor series
% expansion for rotation matrix integration.
% Input:
%   wt : rotation angle in radius, 3x1 vector
% output:
%   c_coeff : coefficient for sin term
function c_coeff = rt_exp_cos(wt)
    d_ang = norm(wt);
    c_coeff = (1 - cos(d_ang))/d_ang^2;
end