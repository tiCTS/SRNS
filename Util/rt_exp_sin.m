%% RT_Expansion_Sin
% This function computes the coefficient of sin term of Taylor series
% expansion for rotation matrix integration.
% Input:
%   wt : rotation angle in radius, 3x1 vector
% output:
%   s_coeff : coefficient for sin term
function s_coeff = rt_exp_sin(wt)
    d_ang = norm(wt);
    s_coeff = sin(d_ang)/d_ang;
end