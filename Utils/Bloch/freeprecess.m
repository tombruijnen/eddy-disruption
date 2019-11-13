function [Afp,Bfp]=freeprecess(T,T1,T2,df)
%
%	Function simulates free precession and decay
%	over a time interval T, given relaxation times T1 and T2
%	and off-resonance df.  Times in ms, off-resonance in Hz.

phi = 2*pi*df*T;	% Resonant precession, radians.
E1 = exp(-T/T1);	
E2 = exp(-T/T2);

Afp = [E2 0 0;0 E2 0;0 0 E1]*zrot(phi);
Bfp = [0 0 1-E1]';


