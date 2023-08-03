function p_next = PositionPropagate(p_ECEF_0, v0, v_next, dt)
% Position integration using RK4
pk1 = v0;
pk2 = (v0 + v_next)/2;
pk3 = (v0 + v_next)/2;
pk4 = v_next;
p_next = p_ECEF_0 + dt/6 * (pk1 +2*pk2 + 2*pk3 + pk4);