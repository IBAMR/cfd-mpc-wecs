%% Calculate the pressure field at location (x, z) at time t:

function [pressure] = calculate_pressure(x, t, z)

global wave hydro

if ( strcmp(wave.wave_type,'FIRST_ORDER_STOKES') )
    
    pressure = - 0.0*wave.rho_w * wave.gz * z + wave.rho_w * wave.gz * (wave.H_wave/2) * cosh(wave.wave_no*(wave.depth + z)) .* cos(wave.wave_no * x - wave.omega * t) / cosh(wave.wave_no * wave.depth);

elseif ( strcmp(wave.wave_type,'IRREGULAR') )

    N = wave.N_components;
    pressure = 0.0;
    
    for m = 1 : N

        pressure = pressure + wave.rho_w * wave.gz * wave.Amp(m) * cosh(wave.k_no(m)*(wave.depth + z)) .* cos(wave.k_no(m) * x - wave.Ohm(m) * t + wave.Phase(m)) / cosh(wave.k_no(m) * wave.depth);
     
    end

    pressure = pressure - 0.0*wave.rho_w * wave.gz * z;
end

end