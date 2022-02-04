function [P, Ju, Jv] =  mpc_matrix(A, B, C, F, Np, Nc)

[Cn, Cm] = size(C);
[An, Am] = size(A);
[Bn, Bm] = size(B);
[Fn, Fm] = size(F);

P  = zeros(Np*Cn, Am);
Ju = zeros(Np*Cn, Nc*Bm);
Jv = zeros(Np*Cn, Np*Fm);
h  = zeros(Np*Cn, Cm);


P(1:Cn,:) = C*A; 
h(1:Cn,:) = C;
for kk = 2:Np
    P((kk-1)*Cn+1:kk*Cn, :) = P((kk-2)*Cn+1:(kk-1)*Cn, :)*A;
    h((kk-1)*Cn+1:kk*Cn, :) = P((kk-2)*Cn+1:(kk-1)*Cn, :);
end

v = h*B;
Ju(:,1:Bm) = v;
for kk = 2:Nc
    Ju(:,(kk-1)*Bm+1 : kk*Bm) = [zeros(Cn*(kk-1), Bm); v(1:Cn*(Np-kk+1),:)]; 
end

q = h*F;
Jv(:, 1:Fm) = q;
for kk = 2:Np
    Jv(:,(kk-1)*Fm+1 : kk*Fm) = [zeros(Cn*(kk-1), Fm); q(1:Cn*(Np-kk+1),:)]; 
end