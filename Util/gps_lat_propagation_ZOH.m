function gps_lat_next = gps_lat_propagation_ZOH(vel0, lat0, h0, dt)
    d_lat = vel0(2)/(M_LLF(lat0)+h0);

    gps_lat_next = lat0 + d_lat * dt;