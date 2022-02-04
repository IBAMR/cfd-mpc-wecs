function F_exc = F_excitation(t)

global hydro wave

if ( strcmp(wave.wave_type,'FIRST_ORDER_STOKES') )
    F_exc = Fexc_for_1st_order_waves(t);
elseif ( strcmp(wave.wave_type,'IRREGULAR') )     
    F_exc = Fexc_for_irreg_waves(t);
else
    msg = 'wave_type needs to be one of FIRST_ORDER_STOKES or IRREGULAR';
    error(msg)
end

end