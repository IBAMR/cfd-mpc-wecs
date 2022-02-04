%% Calculating tyhe NLFK force:
%
%  F_diff(t) = - integral{-tf}{tf}( K_d(tau) . eta(t-tau;x_B) dtau )
%

function F_diff = calculate_diffraction_force(time)

global wave hydro

tau = hydro.tau_diff1;
K_diff = hydro.K_diff1;

y = zeros(length(tau),length(time));
eta_exact = zeros(length(tau),length(time));
Fdiff_exact = zeros(1,length(time));
t = zeros(1,length(tau));

for i = 1 : length(time)
   for j = 1 : length(tau)
       t(j) = tau(j);
       eta_exact(j, i) = (wave.H_wave/2)                      ...
                        * cos( wave.wave_no * hydro.x_B       ...
                             - wave.omega * (time(i) - tau(j)) );
       y(j, i) = eta_exact(j, i) * K_diff(1, j);
   end
   Fdiff_exact(i) = trapz(t, y(:, i)) * wave.rho_w * hydro.g; 
end

F_diff = Fdiff_exact';

end