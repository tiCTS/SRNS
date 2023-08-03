function lat_next = gps_lat_propagation_FOH(v0, v1, lat0, h0, M, dt)

lat_next = lat0 + 0.5*(v0(2) + v1(2))*dt/(M+h0);