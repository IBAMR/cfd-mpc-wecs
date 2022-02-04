function [H,Q] = mpc_quadprog_H_matrix(Ju, Np)

Q = zeros(3*Np, 3*Np);

for ii = 1: Np - 1
    Q(3*(ii - 1)+2, 3*(ii - 1)+3) = 1; 
    Q(3*(ii - 1)+3, 3*(ii - 1)+2) = 1;
end
for ii = Np: Np
    Q(3*(ii - 1)+2, 3*(ii - 1)+3) = 1/2; 
    Q(3*(ii - 1)+3, 3*(ii - 1)+2) = 1/2;
end

H = Ju'*Q*Ju;       