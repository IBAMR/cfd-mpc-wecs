function [ S, Amp, Phase ] = JONSWAP(Ohm, Hs, Tp)
%JONSWAP - Calculates the wave spectrum values for a JONSWAP spectrum

wp = 2*pi/Tp;
Gamma = 3.3;
 for x = 1:length(Ohm)
     if Ohm(x)<wp
         Sigma = 0.07;
     else
         Sigma = 0.09;
     end
     A = exp(-((Ohm(x)/wp-1)/(Sigma*sqrt(2)))^2);
     S(x) = 320*Hs^2*Ohm(x)^-5/Tp^4*exp(-1950*Ohm(x)^-4/Tp^4)*Gamma^A;
 end

 % Determine the frequency step from the frequency vector. Note that the
 % highest frequency step is extrapolated.
 domg = zeros( size(Ohm) );
 domg(1:end-1) = diff( Ohm );
 domg(end) = domg(end-1);

 % Determine the amplitudes from the spectral values
 Amp = sqrt( 2 * S .* domg );

 % Random phases
 Phase = rand(1,length(Ohm))*2*pi;
    
end
