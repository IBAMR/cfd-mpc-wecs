%% Load the BEM data obtained from BEMIO.m file, set wave and MPC parameters

global hydro wave          

load('BEM_data.mat');
hydro.heave_dof = 3;         % For heave DOF

%% WEC body parameters:

k_stiffness = hydro.C(hydro.heave_dof,hydro.heave_dof) * (hydro.g*hydro.rho);% rho*g multiplied because C is normalized quantity 
m_inf       = hydro.Ainf(hydro.heave_dof,hydro.heave_dof) * hydro.rho;       % rho multiplied because Ainf is normalized quantity
a           = (mass + m_inf);

%% Viscous force parameters:

hydro.Cd         = 0.000006865;  % Coefficient of drag
hydro.cross_area = pi * r_cyl^2; % Cross sectional area of cylinder
z_dot            = 0.0;

%% MPC matrix:

mpc_start_time = 10 * wave.Tp;             % Time to start MPC        
dt_controller  = 0.05;                     % sampling time (the controller time-step) 

no_of_Tp = 1.0;
prediction_time = no_of_Tp * wave.Tp;      % Prediction horizon

Np = floor(prediction_time/dt_controller); % Prediction horizon points
Nc = floor(prediction_time/dt_controller); % Control horizon points (usually same as Np)

%% Convolution integral parameters used to calculate the excitation force:

wave.tf   = 1.7;                                                                       
wave.df   = wave.omega * wave.tf / wave.wave_no;    % Using, k*df = w*tf to calculate distance between x_A and x_B

t_exc1      = hydro.ex_t;                           % For excitation force
K_exc1(1,:) = hydro.ex_K(hydro.heave_dof, 1, :);    

t_diff      = hydro.sc_t;                           % For diffraction force
K_diff(1,:) = hydro.sc_K(hydro.heave_dof, 1, :);    

hydro.tau1  = -wave.tf : abs(t_exc1(2) - t_exc1(1)) : wave.tf; % tau variable for eta_wave at point x_B
hydro.K_e1  = spline(t_exc1, K_exc1, hydro.tau1);              % K_e not shifted by tf

hydro.tau_diff1  = -wave.tf : abs(t_diff(2) - t_diff(1)) : wave.tf; % tau variable for eta_wave at point x_B
hydro.K_diff1  = spline(t_diff, K_diff, hydro.tau_diff1);           % K_diff not shifted by tf

past_time         = 2 * wave.tf;         % past data of eta_wave at point x_A to be used for Fexc calculation.
dt_sampling       = dt_controller/2;     % sampling time for past data collection
sample_size       = floor(past_time / dt_sampling) + 1; 
sampling_interval = floor(dt_sampling / dt_solver);

hydro.tau_s  = 0 : dt_sampling : 2*wave.tf;                   % tau variable foreta_wave at point x_A.
hydro.K_es  = spline(t_exc1+wave.tf, K_exc1, hydro.tau_s);    % shifting K_exc to the right by tf

%% Wave parameters:

wave.wave_type = 'FIRST_ORDER_STOKES';  % set wave_type to 'FIRST_ORDER_STOKES' or 'IRREGULAR'

if ( strcmp(wave.wave_type,'IRREGULAR') )
    
    if (isfile('irregular_wave.txt'))   % irregular_wave.txt file is generated from CFD code.
        
        disp('Generating irregular waves based on given data in irregular_wave.txt');
        
        irreg_wave_data = load('irregular_wave.txt');
        [N1, M1] = size(irreg_wave_data);
        wave.N_components = N1;
        
        Amp       = irreg_wave_data(:,1);
        wave.Ohm  = irreg_wave_data(:,2);
        wave.k_no = irreg_wave_data(:,3);
        Phase     = irreg_wave_data(:,4);
        
    else
    
        omega_begin = 2 * pi / 5.2358;      
        omega_end   = 2 * pi / 0.3141;

        wave.N_components = 50;               % number of wave components
        wave.Ohm = linspace(omega_begin, omega_end, wave.N_components);  % N_components number of omega values
    
        alpha =  (wave.depth/hydro.g)*wave.Ohm.^2;
        beta  =  alpha.*tanh(alpha).^(-0.5);
        kd    = (alpha + beta.^2.*cosh(beta).^(-2))./(tanh(beta) + beta.*cosh(beta).^(-2));
        wave.k_no = kd/wave.depth;

        seed = 12345;     rng(seed);          % random number generator
    
        wave.wave_spectrum = 'BRETSCHNEIDER'; % set wave_spectrum to 'JONSWAP' or 'BRETSCHNEIDER'
        Hs = wave.H_wave;
        Ts = wave.Tp;
    
        if ( strcmp(wave.wave_spectrum,'JONSWAP') )
            [ Sw, Amp, Phase ] = JONSWAP(wave.Ohm, Hs, Ts);
        elseif ( strcmp(wave.wave_spectrum,'BRETSCHNEIDER') )
            [ Sw, Amp, Phase ] = BRETSCHNEIDER(wave.Ohm, Hs, Ts);
        else
            msg = 'For IRREGULAR wave, wave_spectrum has to be either JONSWAP or BRETSCHNEIDER';
            error(msg)
        end
    end
    
    wave.Amp   = Amp;
    wave.Phase = Phase;
end

%% LFK and NLFK force calculation options:

wave.Fexc_type = 'LINEAR_FK';     % Set Fexc_type to 'LINEAR_FK' or 'NON_LINEAR_FK'
                                  % When using 'NON_LINEAR_FK', set
                                  % AR_start_time to large values as NLFK
                                  % method uses analytical expressions for
                                  % F_exc calculation.
                                  
dX = [ 0.02 ; 0.02 ; 0.02 ]; % dX = [ dx ; dy ; dz ] --> mesh size of static Cartesian grid for NLFK force calculation

%% AR model parameters:
                     
AR_order = 3;                                % order of the AR-model
AR_start_time = 1 * mpc_start_time;          % Use AR-model from this time onwards. 
                                             % Set large value for switching it off 
                                             % and use analytical
                                             % expressions.

AR_past_time         = 2 * wave.Tp;          % past data of eta_wave at point x_A
AR_dt_sampling       = dt_controller/2;      % Sampling time for past data collection 
                                             % to be used in AR predictions
AR_sample_size       = floor(AR_past_time / AR_dt_sampling);
AR_sampling_interval = floor(AR_dt_sampling / dt_solver);

%% Ramp force parameters:
% Ramp function can be used to apply the control force by smoothly 
% increasing its magnitude when MPC starts

t_ramp_start = mpc_start_time;
t_ramp_end   = mpc_start_time + wave.Tp * 0.001; % Ramp function can be enabled 
                                                 % by changing the multiplying factor 0.001

%% Radiation damping variable:

xr_old = zeros(1,hydro.order_max);

%% Constrained MPC:
%  ymin and ymax = [  z  ;  z_dot  ;  control_force/a  ];

ymin = [   0   ;   0   ;  -100/a   ];
ymax = [   0   ;   0   ;   100/a   ];

%% optimization method:
% optimization_methods are 'linear' and 'non-linear'.

optimization_method = 'linear';
    
%% quadprog options:
% For linear optimization
% Refer -> J. Cretel, G. Lightbody, G. Thomas, A. Lewis,
% Maximisation of Energy Capture by a Wave-Energy Point Absorber using Model Predictive Control,
% IFAC Proceedings Volumes 44 (1) (2011) 3714–3721, 18th IFAC World Congress
 
Lambda1 = 2.0;                  % reduces the aggressiveness of the controller
Lambda2 = 0.2;                  % reduces negative PTO power
options_quadprog = optimoptions( @quadprog                            ...
                               , 'Algorithm'           , 'interior-point-convex'...
                               , 'MaxIterations'       , 5000         ...
                               , 'ConstraintTolerance' , 1.0e-6       ...
                               , 'StepTolerance'       , 1.0e-6       ...
                               , 'OptimalityTolerance' , 1.0e-6    );
                  
%% fmincon options:
% non-linear optimization 
% Refer -> Paolino Tona, Hoaï-Nam Nguyen, Guillaume Sabiron, Yann Creff. 
% An Eﬀiciency-Aware Model Predictive Control Strategy for a Heaving Buoy Wave Energy Converter. 
% 11th European Wave and Tidal Energy Conference - EWTEC 2015, Sep 2015, Nantes, France. hal-01443855

eta_0  = 1.0;                  % reduces negative PTO power
lambda = 7;                    % reduces the aggressiveness of the controller
options_fmincon = optimoptions( @fmincon                              ...
                              , 'SpecifyObjectiveGradient'  ,  true   ...
                              , 'Algorithm'                 , 'sqp'   ...
                              , 'MaxIterations'       , 5000          ...
                              , 'ConstraintTolerance' , 1.0e-6        ...
                              , 'StepTolerance'       , 1.0e-6        ...
                              , 'OptimalityTolerance' , 1.0e-6    );
                          
deltaU = zeros(Np, 1);        % Initial condition 
