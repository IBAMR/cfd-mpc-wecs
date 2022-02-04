%% To simulate a vertical cylindrical point absorber with MPC for heave DOF:

clc;   clear;   close all;

%% Cummins equation:   
% 
%   z_dot_dot(t) = ( 1/(M + A_inf) ) * ( F_h + F_r + F_exc + F_v + F_control )
%
% For heave dof: z(t)         -> dispacement
%                z_dot(t)     -> velocity
%                z_dot_dot(t) -> acceleration 
%
%   X_dot(t) = A . X(t) + B . ( u(t) + v(t) )  -->  Dynamics
%   z(t)     = C . X(t)                        -->  Measurement (position and velocity)
% 
%   where,
%   X = [ z(t)  z_dot(t) x1  x2  x3 ... xr ]' 
%   where, x1,...,xr -> from state-space approx. of the convolution integral
%                    for radiation force
%          r    -> order_max       -->  Order of the state-space system for calculating radiation force  
%          v(t) -> F_exc(t)/a      -->  Excitation force per unit mass
%          u(t) -> F_control(t)/a  -->  Control action per unit mass
%   
%          |                                        |
%          |        0           1          0        |
%          |                                        |
%    A  =  |  -(k_stiff/a)   2*beta/a   -Cr(1xn)/a  | where,   
%          |                                        |    a       -> (m + m_inf)
%          |        0         Br(nx1)    Ar(nxn)    |    k_stiff -> hydrostatic stiffness
%          |                                        |    beta -> -(1/2) * Cd * rho * A_cross * z_dot(0)
%                                                        Ar   -> ( n x n )
%                                                        Br   -> ( n x 1 )
%                                                        Cr   -> ( 1 x n )
%
%   B = [  0     1     0  ]'
%
%   C = |  1     0     0  |     
%       |  0     1     0  |
%

%% Solver parameters:

dt_solver = 0.002;  % Solver time-step
time_end  = 40;     % End time for the simulation

%% Set wave parameters and call script to load .mat file and set MPC parameters:

addpath([cd,'/MPC_matlab_code']);           % Setting path to MPC code directory

global hydro wave mpc_interface

wave.H_wave     = 0.5;                      % Wave height
wave.Tp         = 1.5652;                   % Time period of the wave
wave.omega      = 2 * pi / wave.Tp;         % Angular frequency of the wave
wave.wavelength = 3.8144;                   % Wavelength of the wave
wave.wave_no    = 2 * pi / wave.wavelength; % Wavenumber of the wave
wave.depth      = 2.0;
wave.rho_w      = 1025;                     % Density of water
wave.gz         = -9.81;

r_cyl       = 0.25;                         % Radius of cylinder    
L_cyl       = 0.8;                          % Length  of cylinder
rho_s       = wave.rho_w / 2.0;             % Density of half submerged cylinder 

mass = pi * r_cyl^2 * L_cyl * rho_s;        % Mass of the cylinder

z_initial = 0.0;                            % Initial z-location of COG of cylinder 

load_mpc_parameters;                        

%% Calculate body location x_B and probe point location x_A:

hydro.x_B = wave.wavelength + 5 * r_cyl;    % x-location of cylinder COG
hydro.x_A = hydro.x_B - wave.df;            % probe point x_A location

%% Create static Cartesian grid for NLFK calculations around the cylinder:

X_0 = [ hydro.x_B  ;  2  ; wave.depth +  z_initial ];  % X_0 = [ X_CG  ;  Y_CG  ;  Z_CG ] --> COG of cylinder

create_mesh;

%% Solver parameters to be used:

mpc_interface.m_plus_Ainf      = a;              % ( m + A_inf ) = ( m + m_inf )
mpc_interface.dt_controller    = dt_controller;
mpc_interface.dt_cfd           = dt_solver;      
mpc_interface.mpc_start_time   = mpc_start_time;
mpc_interface.initial_position = z_initial;
mpc_interface.current_position = z_initial;

mpc_interface.next_tuning_time = mpc_start_time;
mpc_interface.old_tuning_time  = 0.0;

mpc_interface.dt_cfd_current   = 0.0;

mpc_interface.velocity         = 0.0;
mpc_interface.Fexc_current     = 0.0;
mpc_interface.F_control        = 0.0;
mpc_interface.F_control_old    = 0.0;
mpc_interface.deltaU           = 0.0;
mpc_interface.eta_A_current    = 0.0;

mpc_interface.dt_sampling       = dt_sampling;
mpc_interface.sample_size       = sample_size;
mpc_interface.sampling_interval = sampling_interval;

%% Time loop:

mpc_interface.t_past = zeros(mpc_interface.sample_size, 1);
mpc_interface.eta_A  = zeros(mpc_interface.sample_size, 1);

AR_t_past  = zeros(AR_sample_size, 1);
AR_eta_A   = zeros(AR_sample_size, 1);
AR_z       = zeros(AR_sample_size, 1);

t_past     = zeros(mpc_interface.sample_size, 1);
eta_A_past = zeros(mpc_interface.sample_size, 1);

nsteps     = floor(time_end / dt_solver);
t_discrete = zeros(nsteps, 1);
y_discrete = zeros(nsteps,hydro.order_max+2);
f_discrete = zeros(nsteps, 1);
Fv_linearized = zeros(nsteps, 1);
Fpto_discrete = zeros(nsteps, 1);
eta_exact  = zeros(nsteps, 1); 
ud         = 0;
        
y_discrete(1,1) = z_initial;

time_current  = 0;
iteration_num = 0;
k = 0;

while (time_current <= time_end)

    k = k + 1;
    time_current
     
    dt_cfd = dt_solver;
    z      = mpc_interface.current_position;
    z_dot  = mpc_interface.velocity;
    
    calculate_mpc_matrices;
    get_radiation_damping_xr;
    
    if ( time_current >= mpc_interface.next_tuning_time ) % Initially, next_tuning_time = mpc_start_time 

        F_control = mpc_interface.F_control;
        mpc_interface.F_control_old = mpc_interface.F_control;
        
        t_past = mpc_interface.t_past;
        eta_A_past = mpc_interface.eta_A;
        
        t_past(end) = time_current;
        eta_A_past(end) = mpc_interface.eta_A_current;
         
        AR_eta_A_past = AR_eta_A;
        AR_z_past     = AR_z;
        
        AR_t_past(end) = time_current;
        AR_eta_A_past(end) = mpc_interface.eta_A_current;
        
        get_control_force;

        mpc_interface.deltaU = deltaU_first;

        mpc_interface.old_tuning_time  = mpc_interface.next_tuning_time;
        mpc_interface.next_tuning_time = mpc_interface.old_tuning_time + mpc_interface.dt_controller;

        fprintf('Next tuning time = %f s', mpc_interface.next_tuning_time);
          
    end
    
    mpc_interface.F_control = mpc_interface.F_control_old                        ...
                            +  mpc_interface.m_plus_Ainf * mpc_interface.deltaU  ...
                            * (time_current + dt_solver - mpc_interface.old_tuning_time)/mpc_interface.dt_controller; % Calculating F_control at current simulation time
    
    control_force = mpc_interface.F_control;
    
    if ( strcmp(wave.Fexc_type,'LINEAR_FK') )
        
        Fexc = F_excitation(time_current + dt_solver);  % Calculating F_exc as LFK at current simulation time
        
    elseif (strcmp(wave.Fexc_type,'NON_LINEAR_FK'))
        
        if ( k <= 10 )
            z_new = z + dt_cfd * z_dot; % Approx. new z-location of the cylinder
        else
            z_new = (1/3)*(2 * dt_cfd * z_dot + 4*y_discrete(k,1) - y_discrete(k-1,1)); % Approx. new z-location of the cylinder
        end

        X_0   = [ hydro.x_B  ;  2  ; wave.depth + z_new ];                    % New COG of the cylinder 
        [psi] = get_level_set_for_cylinder(var(:,1:3), X_0, r_cyl, L_cyl);    % Level-set for cylinder
        [phi] = get_level_set_for_wave(var(:,1:3), time_current + dt_solver); % Level-set for waves

        var(:, 4) = psi';
        var(:, 5) = phi;
        
        F_diff(k+1,1) = calculate_diffraction_force(time_current + dt_solver);   % Calculating difffracion force at current simulation time
        F_fk(k+1,:)   = calculate_NLFK_force(var, time_current + dt_solver, dX); % Calculating NLFK force at current simulation time
        
        Fexc = F_diff(k+1,1) + F_fk(k+1,3); % Calculating F_exc as (NLFK+Diff) at current simulation time
        
    end
    
    f_discrete(k+1,1) = Fexc; 
    
    beta1 = (-1/2) * hydro.Cd * hydro.rho * hydro.cross_area;
    Fv0(k+1,1) = - beta1 * abs(z_dot) * z_dot;  % From linearized part of Fv
             
    y_old = [  y_discrete(k,1:2)  ,  xr_current'  ];
    
    y_new = y_old' + dt_solver * ( Ac * y_old'    ...
                                 + Bc * (f_discrete(k+1,1) + control_force + Fv0(k+1,1))/a ); % Advancing the dynamics
       
    y_discrete(k+1,1:hydro.order_max+2) =  y_new';
    t_discrete(k+1,1)                   = time_current + dt_solver;
    Fpto_discrete(k+1,1)                = control_force;
    
    if ( strcmp(wave.wave_type,'FIRST_ORDER_STOKES') )
        eta_exact(k+1, 1) = (wave.H_wave/2)                             ...
                           * cos( wave.wave_no * hydro.x_A              ...
                                - wave.omega * (time_current+dt_solver) ); % Wave elevation for 1st order waves
    
    elseif ( strcmp(wave.wave_type,'IRREGULAR') )
        eta = 0.0;
        for m = 1:wave.N_components      
            eta = eta + wave.Amp(m)*cos(wave.k_no(m)*hydro.x_A - wave.Ohm(m)*(time_current+dt_solver) + wave.Phase(m));
        end
        eta_exact(k+1, 1) = eta; % Wave elevation for irregular waves
    end
    
    mpc_interface.current_position = y_discrete(k+1,1);           
    mpc_interface.velocity         = y_discrete(k+1,2);            
    mpc_interface.Fexc_current     = f_discrete(k+1,1);
    mpc_interface.eta_A_current    = eta_exact(k+1,1);
    
    Fv_linearized(k+1,1) = 2 * beta1 * abs(z_dot) * y_discrete(k+1,2)...
                         - beta1 * abs(z_dot) * z_dot;         
   
    if ( (mpc_interface.t_past(end) + mpc_interface.dt_sampling) <= (time_current+dt_solver) ) % Collecting past wave data
                    
        fprintf('mpc t_past = %f s', mpc_interface.t_past(end) + mpc_interface.dt_sampling);
        
        mpc_interface.t_past(1:end-1) = mpc_interface.t_past(2:end);
        mpc_interface.eta_A(1:end-1)  = mpc_interface.eta_A(2:end);
        
        mpc_interface.t_past(end) = t_discrete(k+1,1);
        mpc_interface.eta_A(end)  = mpc_interface.eta_A_current;
                   
    end
    
    if ( (AR_t_past(end) + AR_dt_sampling) <= (time_current+dt_solver) ) % Collecting past z-location data of the cylinder
        
        fprintf('AR t_past = %f s', AR_t_past(end) + AR_dt_sampling);
        
        AR_t_past(1:end-1) = AR_t_past(2:end);
        AR_eta_A(1:end-1)  = AR_eta_A(2:end);
        AR_z(1:end-1)      = AR_z(2:end);
        
        AR_t_past(end) = t_discrete(k+1,1);
        AR_eta_A(end)  = mpc_interface.eta_A_current;
        AR_z(end)      = mpc_interface.current_position;
                   
    end

    iteration_num = iteration_num + 1;       % Iteration counter
    time_current = time_current + dt_solver; % Advance simulation time
    
end

%% Save the final result:

save('BEM_results.mat','t_discrete','y_discrete','Fpto_discrete','f_discrete','eta_exact');
        
%% Plot the results:

figure(1);
plot(t_discrete, y_discrete(:,1));       % z vs t plot
ylabel ('Displacement (m)');
xlabel ('Time (s)');
grid on;

figure(2);
yyaxis left;
plot(t_discrete, y_discrete(:,2), 'b-'); % z_dot vs t plot
ylabel('Velocity (m/s)') 
grid on;

yyaxis right;
plot(t_discrete, f_discrete, 'r-');      % F_exc vs t plot
ylabel('Excitation force (N)','Color','r') 
grid on;
xlabel ('Time (s)');

figure (3);
plot(t_discrete, Fpto_discrete, 'r-');   % F_PTO vs t plot
ylabel ('Control force (N)');
xlabel ('Time (s)');
grid on;

avg_start = (30);                        % Calculate avg. power from t = 30 s and onwards  
avg_start_step = floor(avg_start/dt_solver);
avg_pto_power = mean(-Fpto_discrete(avg_start_step:end).*y_discrete(avg_start_step:end,2)) % P_avg value

figure(4);
hold on;
plot(t_discrete, -Fpto_discrete.*y_discrete(:,2));                 % P vs t plot
plot(t_discrete, avg_pto_power*ones(1,length(t_discrete)),'--k' ); % P_avg vs t plot
hold off;
ylabel ('PTO power (W)');
xlabel ('Time (s)');
grid on;
