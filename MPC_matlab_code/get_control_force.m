%% To get control force:

global hydro wave

time_horizon_shift = time_current + (1:Np)' * dt_controller - dt_controller;  % t = t_current to (t_current+Th-dt_controller)
time_horizon       = time_current + (1:Np)' * dt_controller;                  % t = (t_current+dt_controller) to (t_current+Th)
time1 = [ time_horizon_shift ; time_horizon_shift(end) + dt_controller ];     % t = t_current to (t_current + Th)
Nt = length(time1);            

if ( time_current < AR_start_time ) 
            
    disp('Using analytical expression for F_excitation calculation');
    
    if (strcmp(wave.Fexc_type,'LINEAR_FK'))
    
        F_exc1 = F_excitation(time_horizon_shift);                      % F_exc for t = t_current to (t_current+Th-dt_controller)
        F2     = F_excitation(time_horizon_shift(end) + dt_controller); 
        F_exc2 = [  F_exc1(2:end)  ;  F2  ];                            % F_exc for t = (t_current+dt_controller) to (t_current+Th)
    
    elseif (strcmp(wave.Fexc_type,'NON_LINEAR_FK'))
        
        z_pred = AR_prediction(AR_z_past, AR_t_past, Np, AR_order, dt_controller); % Using AR prediction for z over the time horizon
        F_diff1 = calculate_diffraction_force(time1);  % Diffracation force calculation over time horizon
        
        F_fk1 = zeros(Nt,3);
        for m = 1 : Nt
            
            X_0(3) = wave.depth + z_pred(m);
            [psi]  = get_level_set_for_cylinder(var, X_0, r_cyl, L_cyl); % Levet-set for cylinder
            [phi]  = get_level_set_for_wave(var, time1(m));              % Level-set for wave
            var(:, 4) = psi';
            var(:, 5) = phi;
            
            F_fk1(m,:) = calculate_NLFK_force(var, time1(m), dX); % Calculating NLFK force
        
        end
               
        F_exc1(:,1) = F_diff1(1:end-1,1) + F_fk1(1:end-1,3); % F_exc for t = t_current to (t_current+Th-dt_controller)
        F_exc2(:,1) = F_diff1(2:end,1) + F_fk1(2:end,3);     % F_exc for t = (t_current+dt_controller) to (t_current+Th)
        
    end

    deltaF = F_exc2 - F_exc1;
    
else 
    
    disp('Using AR model for F_excitation calculation');
            
    y_predicted = AR_prediction(AR_eta_A_past, AR_t_past, Np, AR_order, dt_controller); % Using AR prediction to get 
                                                                                        % eta_wave at point x_A in time horizon
     
    tau_s = hydro.tau_s;
    K_e = hydro.K_es;       % Using shifted IRF for wave excitation
    
    TIME  = [ time_current  ;  time_horizon ];                      % TIME = t to (t+Tp).
    t1    = time_current : dt_sampling : time_horizon(end);                
    y_p   = spline([time_current ; time_horizon], y_predicted, t1); % Adjusting data with dt_sampling interval
    eta_A = [ eta_A_past    ;   y_p(2:end)' ];                      % eta_A -> (t-2*tf) to (t+Th).
    
    t_start_idx = length(eta_A_past);   
        
    for i = 1 : length(TIME) 
        for j = 1 : length(tau_s)
            if ( tau_s(j) <= 2*wave.tf )
                t2(j) = tau_s(j);
                eta_A_idx = (i - 1 + t_start_idx) - (j - 1);
                eta1(j, i) = eta_A(eta_A_idx);
                y(j, i) = eta1(j, i) * K_e(j);
            end
        end
        F_pred(i) = trapz(t2, y(:, i)) * hydro.rho * hydro.g; % Multiplying my (rho_w*g) because K_e is normalized quantity
    end
   
    F_exc1 = F_pred(1:end-1); % F_exc for t = t_current to (t_current+Th-dt_controller)
    F_exc2 = F_pred(2:end);   % F_exc for t = (t_current+dt_controller) to (t_current+Th)

    deltaF = F_exc2' - F_exc1';
    
end

F_v = - (-1/2) * hydro.Cd * hydro.rho * hydro.cross_area * abs(z_dot) * z_dot; % From linearized part

Fexc_current = F_exc1(1) + F_v;

x = [  z  ;  z_dot  ;   xr_current  ;  F_control/a  ;  Fexc_current/a  ];

G = mpc_quadprog_G_matrix(P, Ju, Jv, Q, x, deltaF/a);

% Constrained optimization: If no constraints are provided then
% unconstrained optimization is done.
Aineq = [];
Bineq = [];
[Aineq, Bineq] = constraint(P, Ju, Jv, x, deltaF/a, ymin, ymax);

if ( strcmp(optimization_method,'linear') )
    
    disp('Using quadprog for optimization');
    H1 = Lambda1 * eye(size(H));                                     % reduces aggressiveness of controller
    
    L = tril( ones(length(deltaF)) );
    H2 = Lambda2 * (L' * L);                                         % reduces negative absorbed power
    G2 = Lambda2 * 2 * L' * F_control/a * ones(length(deltaF),1);    % reduces negative absorbed power

    deltaU = quadprog(H + 2*H1 + 2*H2, G + G2, Aineq, Bineq, [], [], [], [], [], options_quadprog); % Calculating optimal control sequence

elseif ( strcmp(optimization_method,'non-linear') )
    
    disp('Using efficiency based optimization');
    Px_Jv = P*x + Jv*(deltaF/a);

    deltaU = fmincon(@(du) get_obj_fun_and_grad_fun(du, x, (F_exc2 + F_v)/a, eta_0, Ju, Px_Jv, lambda) ,deltaU , Aineq, Bineq, [], [], [], [], [], options_fmincon); % Calculating optimal control sequence

else
    msg = 'Atleast one optimization method needs to be specified';
    error(msg)
end

deltaU_first = deltaU(1,1); % First element of the optimal control sequence
