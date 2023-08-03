function gps_long_next = gps_long_propagation_ZOH( vel0, long0, lat1, h0, dt)
    d_long = vel0(1) /  ((N_LLF(lat1) + h0)*cos(lat1));

    gps_long_next = long0 + d_long*dt;