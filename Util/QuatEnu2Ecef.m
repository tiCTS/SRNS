function q_LU2E = QuatEnu2Ecef(lat_deg, long_deg)
% generate quaternion from local enu to ecef frame.
% input lat and long are in degree
R_LU2E = [ -sind(long_deg) -cosd(long_deg)*sind(lat_deg) cosd(long_deg)*cosd(lat_deg);
                             cosd(long_deg) -sind(long_deg)*sind(lat_deg) sind(long_deg)*cosd(lat_deg);
                             0 cosd(lat_deg) sind(lat_deg)];


q_LU2E_0 = 0.5*sqrt(1+R_LU2E(1,1)+R_LU2E(2,2)+R_LU2E(3,3));
q_LU2E_1 = 0.5*sqrt(1+R_LU2E(1,1)-R_LU2E(2,2)-R_LU2E(3,3)) * sign(R_LU2E(3,2) - R_LU2E(2,3));
q_LU2E_2 = 0.5*sqrt(1-R_LU2E(1,1)+R_LU2E(2,2)-R_LU2E(3,3)) * sign(R_LU2E(1,3) - R_LU2E(3,1));
q_LU2E_3 = 0.5*sqrt(1-R_LU2E(1,1)-R_LU2E(2,2)+R_LU2E(3,3)) * sign(R_LU2E(2,1) - R_LU2E(1,2));

q_LU2E = quaternion([q_LU2E_0 q_LU2E_1 q_LU2E_2 q_LU2E_3]);