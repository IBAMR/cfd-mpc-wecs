% Calculate the mpc matrices:

[Ar, Br, Cr] = construct_radiation_SS();

%% Constructing the continuous time state-space system:

global Ac Bc Cc

Ac = zeros(hydro.order_max + 2, hydro.order_max + 2);
Bc = zeros(hydro.order_max + 2, 1);
Cc = zeros(2, hydro.order_max + 2);

Ac(1,2) = 1.0;
Ac(2,1) = -k_stiffness/a; 
Ac(2,2) = ( 2 * (-1/2) * hydro.Cd * hydro.cross_area * hydro.rho * abs(z_dot) )/a;
Ac(2,3:hydro.order_max+2) = -(Cr/a) * hydro.rho;   % rho multiplied because Cr is normalized quantity   
Ac(3:hydro.order_max+2,2) = Br;
Ac(3:hydro.order_max+2,3:hydro.order_max+2) = Ar;        
 
Bc(2,1) = 1.0;

Cc(1,1) = 1;
Cc(2,2) = 1; 

%% Constructing the discrete time state-space system using FOH:
%
%   xd[k+1] = Phi.xd[k] + Gamma.(u[k] + v[k]) + Lambda.(d(u)[k+1] + d(v)[k+1])
%
%   yd[k] = Cc.xd[k]
%
%       where, u = F_PTO / (m + m_inf)
%              v = F_exc / (m + m_inf)
%                

[Phi, Gamma, Lambda] = discrete_time_matrix(Ac, Bc, dt_controller);

%% Augmented matrix form:
%
%   x[k+1] = A.x[k] + B.du[k+1] + F.dv[k+1]
%   y[k] = C.x[k]
%
%   where,  
%       |    |
%       | xd |
%   x = | ud |      with, ud = control input
%       | vd |            vd = excitation force
%       |    |
%
%       |    |
%   y = | yd |
%       | ud |
%       |    |

[A, B, F, C] = augmented_matrix(Phi, Gamma, Lambda);

%% MPC matrix:
[P, Ju, Jv] = mpc_matrix(A, B, C, F, Np, Nc);

[H, Q] = mpc_quadprog_H_matrix(Ju, Np); 
H = (H + H')/2;