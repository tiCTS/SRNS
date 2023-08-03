function long_next = gps_long_propagation_FOH(v0, v1, lat1, long0, h0, N, dt )

long_next = long0 + 0.5*(v0(1) + v1(1))*dt/((N+h0)*cos(lat1));