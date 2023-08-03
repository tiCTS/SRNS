function gps_lat_next = gps_lat_propagation_RK4(vel0, vel1, lat0, h0, dt)
    lat_k1 = vel0(2)/(M_LLF(lat0)+h0);
    lat_k2 = 0.5*(vel0(2) + vel1(2)) /(M_LLF(lat0 + 0.5*lat_k1*dt) + h0);
    lat_k3 = 0.5*(vel0(2) + vel1(2)) / (M_LLF(lat0 + 0.5*lat_k2*dt) + h0);
    lat_k4 = vel1(2) / (M_LLF(lat0 + lat_k3*dt) + h0);

    gps_lat_next = lat0 + dt/6 * (lat_k1 + 2*lat_k2 + 2*lat_k2 + lat_k4);
end

