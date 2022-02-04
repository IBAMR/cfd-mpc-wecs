%% Constructing the state-space matrices for solving the radiation force term convolution integral:

function [Ar, Br, Cr] = construct_radiation_SS()
global hydro

Ar = zeros(hydro.order_max,hydro.order_max);
Br = zeros(hydro.order_max,1);
Cr = zeros(1,hydro.order_max);

for i = 1 : hydro.order_max
    for j = 1 : hydro.order_max
        
        Ar(i,j) = hydro.ss_A(hydro.heave_dof,hydro.heave_dof,i,j);
        
    end
    Br(i,1) = hydro.ss_B(hydro.heave_dof,hydro.heave_dof,i);        
    Cr(1,i) = hydro.ss_C(hydro.heave_dof,hydro.heave_dof,1,i);
end