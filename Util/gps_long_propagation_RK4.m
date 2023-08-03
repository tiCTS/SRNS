function gps_long_next = gps_long_propagation_RK4(vel0, vel1, long0, lat1, h0, dt)
    long_k1 = vel0(1) / ((N_LLF(lat1) + h0)*cos(lat1));
    long_k2 = 0.5*( vel0(1) + vel1(1)) / ((N_LLF(lat1) + h0)*cos(lat1));
    long_k3 = 0.5*( vel0(1) + vel1(1)) / ((N_LLF(lat1) + h0)*cos(lat1));
    long_k4 = vel1(1) / ((N_LLF(lat1) + h0)*cos(lat1));

    gps_long_next = long0 + dt/6 * (long_k1 + 2*long_k2 + 2*long_k3 + long_k4);
end

