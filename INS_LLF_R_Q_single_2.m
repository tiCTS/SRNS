% Sounding Rocket Navigation Systems (SRNS)
% INS in Local Level Frame (LLF)
% Using Rotation Matrix & Qiaternion 1st-order hold
% Single precision version
%% Initialization settings

% Q or R or AHRS attitude directly
USE_Q = false;
USE_R = true;
USE_AHRS = false;

SRPL_DATA = false;
GPS_INIT_DATA = true;
GPS_UPDATE = false;

% Update selection%
% Use Kalman Filter
USE_KF = true;
% USE_KF = false;


ZERO_V_UPDATE = true;
% ZERO_V_UPDATE = false;

% AHRS_UPDATE = true;
AHRS_UPDATE = false;

%% parameters
g = single(9.8106);  % defined by AHRS-10P, m/s^2.
a_Earth = single(6378137); % WGS84 Earth semi-major axis (m).
e_sq_Earth = single(6.694379990e-3); % WGS84 Earth eccentricity squared.
w_e =  single(7.292115e-5);  % rad/s

% imu params
sigma_GB = single( 1.932e-07);%1.932e-04; %rad/s
sigma_AB = single(1.954e-04); %1.954e-04;  %m/s^2
TBA = single(1); %sec
TBG_X = single(10);
TBG_Y = single(10);
TBG_Z = single(5);
TBGinv = [1/TBG_X 0 0; 0 1/TBG_Y 0; 0 0 1/TBG_Z];
TBGM = [TBG_X 0 0; 0 TBG_Y 0; 0 0 TBG_Z];

% NCKU SRPL GPS
SRPL_lat = single(22.93767);
SRPL_long = single(120.27274);
SRPL_alt = single(20);

% GPS position & velocity initialization
if GPS_INIT_DATA == true  % already radius
    gps_lat_init = single( gps_init.Lat_init);
    gps_long_init = single( gps_init.Long_init);
    gps_alt_init = single( gps_init.Alt_init);
end
gps_pos(:,1) = [gps_lat_init; gps_long_init; gps_alt_init];
enu_pos(1,:) =single( zeros(1,3));
vel(:,1) = single( [0;0;0] );
vel_RK4(:,1) = single([0;0;0]);   % RK4 comparision

% attitude initialization
q_b2LU = single( quaternion(imu.q0(1), imu.q1(1), imu.q2(1), -imu.q3(1)));
R_b2LU_t1 = Quat2Rot(q_b2LU);
if USE_Q || USE_AHRS
    eulZYX(1,:) = rad2deg(quat2eul(q_b2LU,'ZYX'));
end
if USE_R
    eulZYX(1,:) = rad2deg(rotm2eul(R_b2LU_t1, 'ZYX'));
end
eul_AHRS(1,:) = rad2deg(quat2eul(quaternion(imu.q0(1), imu.q1(1), imu.q2(1), imu.q3(1))));

% imu measurement
w_t1 = single([ imu.GyroX(1); imu.GyroY(1); imu.GyroZ(1)] );
w_t1_LU = QuatVectorTrans(q_b2LU, w_t1);
a_t1 = single( [ imu.AccX(1); imu.AccY(1); imu.AccZ(1)] );
mag_t1 = single([ imu.MagX(1); imu.MagY(1); imu.MagZ(1)]);
if USE_Q || USE_AHRS
    a_t1_LU = QuatVectorTrans(q_b2LU,a_t1);
end
if USE_R
    a_t1_LU = R_b2LU_t1*a_t1;
end
% specific force (SPF)
g_LU(:,1) = GravityFromLLA(gps_pos(1,1),gps_pos(3,1));

% bias
dw(:,1) = single(zeros(3,1));
da(:,1) = single(zeros(3,1));

% KF initialization
del_x = single(zeros(12,1));
% del_x_att = single(zeros(6,1));  % conatains error angle and gyro bias, for attitude estimation.
P = single(1e-3*eye(12));
% P_att = signle()
Q = single(1e-4);
R_zv = single([1e-6*eye(3) zeros(3,4); ...
                          zeros(3,3) 1*eye(3) zeros(3,1); ...
                          zeros(1,6) 1e-2]);

%% Propagate
for idx = 2:size(imu,1)
% for idx = 2:3000
    % claculate params
    M = (a_Earth*(1-e_sq_Earth))/(1 - e_sq_Earth*sin(gps_pos(1,idx-1))^2)^(3/2);
    N = a_Earth/(1 - e_sq_Earth*sin(gps_pos(1,idx-1))^2)^(1/2);

     % get new imu measurement
     w_t2 = single( [ imu.GyroX(idx); imu.GyroY(idx); imu.GyroZ(idx)])- dw(:,idx-1);
%      w_t2_LU = QuatVectorTrans(q_b2LU(idx-1), w_t2);
    w_t2_LU = R_b2LU_t1 * w_t1;
     a_t2 = single([ imu.AccX(idx); imu.AccY(idx); imu.AccZ(idx)]);
    % get dt
     dt = single(imu.SRNS_time(idx) - imu.SRNS_time(idx-1));
     % Attitude propagation
      wE_l = fun_wE_LU(gps_pos(1,idx-1));
      w_el_l = fun_w_el_LU(vel(:,idx-1), M, N, gps_pos(1,idx-1), gps_pos(3,idx-1));
      if USE_Q
          q_b2LU(idx) = Quaternion1stInt(w_t1_LU, w_t2_LU, q_b2LU(idx-1), dt);
          R_b2LU_t2 = Quat2Rot(q_b2LU(idx));
%           [qq0,qq1,qq2,qq3] = parts(q_b2LU(idx));  % for debug
          eulZYX(idx,:) = rad2deg(quat2eul(q_b2LU(idx),'ZYX'));
      end
      if USE_R
          R_b2LU_t2 = rt_integration_LLF(R_b2LU_t1, w_t1_LU, wE_l, w_el_l, dt);
          eulZYX(idx,:) = rad2deg(rotm2eul(R_b2LU_t2, 'ZYX'));
      end
      if USE_AHRS
          q_b2LU(idx) = single( quaternion(imu.q0(idx), imu.q1(idx), imu.q2(idx), imu.q3(idx)));
          R_b2LU_t2 = Quat2Rot(q_b2LU(idx));
          eulZYX(idx,:) = rad2deg(quat2eul(q_b2LU(idx),'ZYX'));
      end
       eul_AHRS(idx,:) = rad2deg(quat2eul(quaternion(imu.q0(idx), imu.q1(idx), imu.q2(idx), imu.q3(idx))));

       % transform acc to LLF
       if USE_Q || USE_AHRS
           a_t2_LU = QuatVectorTrans(q_b2LU(idx), a_t2);
       end
       if USE_R
            a_t2_LU = R_b2LU_t2*a_t2;
       end

        % Velocity propagation
        g_LU(:,idx) = GravityFromLLA(gps_pos(1,idx-1), gps_pos(3,idx-1));
        %RK4
%         vel_RK4(:,idx) = VelocityPropagate_LLF_RK4(vel(:,idx-1), a_t1_LU, a_t2_LU, gps_pos(1,idx-1), gps_pos(3,idx-1), M, N, g_LU(:,idx), dt);
        % ZOH
        vel(:,idx) = VelocityPropagate_LLF_ZOH(vel(:,idx-1), a_t1_LU, gps_pos(1,idx-1), gps_pos(3,idx-1), M, N, g_LU(:,idx), dt);

        % gps position propagation
%         lat_next = gps_lat_propagation_RK4(vel(:,idx-1), vel(:,idx), gps_pos(1,idx-1),gps_pos(3,idx-1),dt);
%         lat_next = gps_lat_propagation_ZOH(vel(:,idx-1), gps_pos(1,idx-1),gps_pos(3,idx-1),dt);
            lat_next = gps_lat_propagation_FOH(vel(:,idx-1), vel(:,idx), gps_pos(1,idx-1), gps_pos(3,idx-1), M, dt);
%         long_next = gps_long_propagation_RK4(vel(:,idx-1), vel(:,idx), gps_pos(2,idx-1), lat_next, gps_pos(3,idx-1), dt);
%         long_next = gps_long_propagation_ZOH(vel(:,idx-1), gps_pos(2,idx-1), lat_next, gps_pos(3,idx-1), dt );    
            long_next = gps_long_propagation_FOH(vel(:,idx-1), vel(:,idx), lat_next, gps_pos(2,idx-1),gps_pos(3,idx-1), N, dt);
%         alt_next = gps_alt_propagation_RK4(vel(:,idx-1), vel(:,idx), gps_pos(3,idx-1), dt);
        alt_next = gps_alt_propagation_ZOH(vel(:,idx-1), gps_pos(3,idx-1), dt);
%             alt_next = gps_alt_propagation_FOH(vel(:,idx-1), vel(:,idx), gps_pos(3,idx-1), dt);

            gps_pos(:,idx) = [lat_next; long_next; alt_next];

        % bias
       dw(:,idx) = dw(:,idx-1) -  TBGinv*dw(:,idx-1)*dt;
        

        if USE_KF
            F_dr_dv = [0 1/(M+gps_pos(3,idx-1)) 0;
                        1/(N+gps_pos(3,idx-1)*cos(gps_pos(1,idx-1))) 0 0;
                        0 0 1];

            F_dv_dang = [0 a_t1_LU(3) -a_t1_LU(2);
                            -a_t1_LU(3) 0 a_t1_LU(1);
                            a_t1_LU(2) -a_t1_LU(1) 0];

            F_dang_dv = [ 0 1/(M+gps_pos(3,idx-1)) 0;  %p
                             -1/(N+gps_pos(3,idx-1))  0 0;  % r
                             -tan(gps_pos(1,idx-1))/(N+gps_pos(3,idx-1))   0 0];  %y

            w_il_l = wE_l + w_el_l;
             F_dang_dang =  -skewMatrix(w_il_l) ;
             F_dang_dw = R_b2LU_t1; 

             F_dw_dw = -TBGinv*eye(3); 
              % F_dw_dw = eye(3);
%                F_dw_dw = zeros(3,3);

              % F structure
    %       dr   dv   dang  dw
    %   dr
    %   dv
    %   dang
    %   dw  

            F = [zeros(3,3) F_dr_dv zeros(3,6);
                    zeros(3,6) F_dv_dang zeros(3,3);
                    zeros(3,3) F_dang_dv F_dang_dang F_dang_dw;
                    zeros(3,9) F_dw_dw];
            G_dw = sqrt( 2*TBGM*(sigma_GB^2)*ones(3,1));
            G = [zeros(9,1); G_dw];

            % State transition matrix
            StateTrans = eye(12) + dt * F ; % 1st-order approx.
            % KF pridiction
            del_x_next = StateTrans * del_x;
            P_next = StateTrans*P*StateTrans' + G*Q*G';

            % Updates include zero V and mag
            if ZERO_V_UPDATE

                mag_meas = [imu.MagX(idx); imu.MagY(idx); imu.MagZ(idx)];
                mag_meas_t1 = [imu.MagX(idx-1); imu.MagY(idx-1); imu.MagZ(idx-1)];

                dR_wt2 = R_b2LU_t2 * R_b2LU_t1' ;
                mag_pre_t2 = dR_wt2*mag_meas_t1;
                del_z_mag = (mag_meas - mag_pre_t2);

                del_z_zv = single([0;0;0]) - vel(:,idx);
                del_z_az = single(eulZYX(1,1) - eulZYX(idx,1));

                del_z = [del_z_zv; del_z_mag; del_z_az];

                               
                H_ZV = single([zeros(3,3) eye(3) zeros(3,6);...
                    zeros(3,6) -skewMatrix(mag_pre_t2)   zeros(3,3); ...
                    zeros(1,8) 1 zeros(1,3)]);

                
                S = (H_ZV*P_next*H_ZV' + R_zv);
                K = P_next*H_ZV'/S;
                del_x_next = del_x_next + K*del_z;
                P_next = (eye(12) - K*H_ZV)*P_next*(eye(12) - K*H_ZV)' + K*R_zv*K';
%                 P_next = (eye(12) - K*H_ZV)*P_next;
                % error injection
                gps_pos(:,idx) = gps_pos(:,idx) + del_x_next(1:3);
                vel(:,idx) = vel(:,idx) + del_x_next(4:6);
%                 del_x_next(9) = -del_x_next(9);
%                 del_ang = [del_x_next(7); del_x_next(8); del_x_next(9)];
%                 dq = quaternion([single(1) 0.5*del_ang']);
%                 dq = normalize(dq);
%                 q_b2LU(idx) = dq*q_b2LU(idx);
%                 [qq0,qq1,qq2,qq3] = parts(q_b2LU(idx));  % for debug
%                 R_b2LU_t2 = Quat2Rot(q_b2LU(idx));
                del_R = (eye(3) + skewMatrix(del_x_next(7:9)));  % equation 5.97 in 2007 book
                R_b2LU_t2 = del_R*R_b2LU_t2;
                dw(:,idx) = del_x_next(10:12);
%                 if idx >= 15000
%                     w_t1 = w_t2;  % JUST CREATE A BREAK POINT HERE FOR DEBUG
%                 end
                % error states reset
                del_x_next = [zeros(9,1); del_x_next(10:12)];
            end

%             if USE_ATT_KF
%                 F_dang_dang = -skewMatrix(R_b2LU*w_t1);
%                  F_dang_dw = R_b2LU_t1;

%                  F = [F_dang_dang F_dang_dw;
%                         zeros(3,6);]

%                  G = [-eye(3) zeros(3,3);
%                          zeros(3,3) eye(3)];


            del_x = del_x_next;
            P = P_next;

        end
            

        w_t1 = w_t2;
        w_t1_LU = w_t2_LU;
        a_t1_LU = a_t2_LU;
        R_b2LU_t1 = R_b2LU_t2;
        
        
        
end
gps_pos(1,:) = rad2deg(gps_pos(1,:));
gps_pos(2,:) = rad2deg(gps_pos(2,:));

%% calculate displacement
for idx = 1: size(imu,1)
% for idx = 1:3000
    enu_pos(idx,:) = lla2enu(double(gps_pos(:,idx)'), double(gps_pos(:,1)'), "ellipsoid");
end
