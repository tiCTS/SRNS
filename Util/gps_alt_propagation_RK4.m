function alt_next = gps_alt_propagation_RK4(vel0, vel1, alt0, dt)
    alt_k1 = vel0(3);
    alt_k2 = 0.5*(vel0(3) + vel1(3));
    alt_k3 = 0.5*(vel0(3) + vel1(3));
    alt_k4 = vel1(3);

    alt_next = alt0 + dt/6*(alt_k1 + 2*alt_k2 + 2*alt_k3 + alt_k4);