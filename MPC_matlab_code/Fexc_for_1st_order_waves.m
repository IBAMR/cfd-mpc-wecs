%% Calculating Fexc using exact expression of wave elevation at x_A using convolution integral:

function F_exc = Fexc_for_1st_order_waves(t)

global hydro wave
  
TIME = t;
tau2 = hydro.tau1;
K_e = hydro.K_e1;

for i = 1 : length(TIME)
   for j = 1 : length(tau2)
       t2(j) = tau2(j);
       eta_exact(j, i) = (wave.H_wave/2)                      ...
                        * cos( wave.wave_no * hydro.x_B       ...
                             - wave.omega * (TIME(i) - tau2(j)) );
       y2(j, i) = eta_exact(j, i) * K_e(1, j);
   end
   Fexc_exact(i) = trapz(t2, y2(:, i)) * 1025 * 9.81; 
end

F_exc = Fexc_exact';

end