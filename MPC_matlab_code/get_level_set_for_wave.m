%% Level-set for tracking a waves:

function [phi] = get_level_set_for_wave(X, t)

global wave 

if ( strcmp(wave.wave_type,'FIRST_ORDER_STOKES') )
        
    eta_exact = (wave.H_wave/2) ...
              * cos( wave.wave_no * X(:,1) - wave.omega * t ); % Wave elevation for 1st order waves
    
elseif ( strcmp(wave.wave_type,'IRREGULAR') )

    N = wave.N_components;
    eta = 0.0;
    
    for m = 1 : N      
        eta = eta + wave.Amp(m)*cos(wave.k_no(m)*X(:,1) - wave.Ohm(m)*t + wave.Phase(m)); % Wave elevation for irregular waves
    end
    eta_exact = eta;

end

distance = -eta_exact + (X(:,3) - wave.depth); % Level-set for tracking waves
phi = distance; 

end