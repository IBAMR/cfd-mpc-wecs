%% Calculating Fexc using solver data of wave elevation at x_A:

function F_exc = Fexc_for_irreg_waves(t)

global hydro wave

TIME = t;
tau2 = hydro.tau1;
K_e = hydro.K_e1;

for i = 1:length(TIME)
    for j = 1 : length(tau2)
        t2(j) = tau2(j);
        eta = 0.0;
        for m = 1:wave.N_components      
           eta = eta + wave.Amp(m)*cos(wave.k_no(m)*hydro.x_B - wave.Ohm(m)*(TIME(i) - tau2(j)) + wave.Phase(m));
        end
        eta_exact(j, i) = eta;
        y2(j, i) = eta_exact(j, i) * K_e(1, j);
    end
    Fexc_exact(i) = trapz(t2, y2(:, i)) * wave.rho_w * hydro.g;
end

F_exc = Fexc_exact';

end