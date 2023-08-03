function w_el_LU = fun_w_el_LU(v, M, N, lat, h)
    w_el_LU = [-v(2)/(M+h);
                        v(1)/(N+h);
                        v(1)*tan(lat)/(N+h)];