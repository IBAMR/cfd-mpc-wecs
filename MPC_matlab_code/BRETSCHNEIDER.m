%% Calculates the wave spectrum values for a BRETSCHNEIDER spectrum

function [ S, Amp2, Phase2 ] = BRETSCHNEIDER(Ohm, Hs, Tp)

wp = 2*pi/Tp;
T1 = 0.772 * Tp;
 for x = 1:length(Ohm)
     S(x) = 173*Hs^2*Ohm(x)^-5/T1^4*exp(-692*Ohm(x)^-4/T1^4);
 end

 % Determine the frequency step from the frequency vector. Note that the
 % highest frequency step is extrapolated.
 domg = zeros( size(Ohm) );
 domg(1:end-1) = diff( Ohm );
 domg(end) = domg(end-1);

 % Determine the amplitudes from the spectral values
 Amp2 = sqrt( 2 * S .* domg );

 % Random phases
 Phase2 = rand(1,length(Ohm))*2*pi;

end
