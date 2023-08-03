function alt_next = gps_alt_propagation_ZOH(vel0, alt0, dt)
    d_alt = vel0(3);

    alt_next = alt0 + d_alt * dt;