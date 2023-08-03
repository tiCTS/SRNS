function alt_next = gps_alt_propagation_FOH(v0, v1, h0, dt)

alt_next = h0 + 0.5*(v0(3) + v1(3))/dt;