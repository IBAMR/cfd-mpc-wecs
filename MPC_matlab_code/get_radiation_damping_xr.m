%% Solve for the radiation damping xr variables:

xr_current = xr_old' + dt_cfd * ( Ar * xr_old' + Br * z_dot); % Advance xr variable by dt_cfd (I.e. dt_solver)

xr_old = xr_current';