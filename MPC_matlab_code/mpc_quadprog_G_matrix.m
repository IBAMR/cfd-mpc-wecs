function G = mpc_quadprog_G_matrix(P, Ju, Jv, Q, x, v)

G = Ju'*Q*(P*x + Jv*v);