% h  : Sampling time
% Ac : Continuous state-space matrix
% Bc : Continuous forcing matrix 

function [Phi, Gamma, Lambda] = discrete_time_matrix(Ac, Bc, h)

[n,m] = size(Ac);
Phi = expm(Ac*h);
Gamma = Ac  \ ((Phi - eye(n,m))*Bc);
Lambda = Ac \ (Gamma - h*Bc);
Lambda = Lambda/h;