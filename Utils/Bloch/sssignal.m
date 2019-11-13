% 
%	function [Msig,Mss] = sssignal(flip,T1,T2,TE,TR,dfreq)
% 
%	Calculate the steady state signal at TE for repeated
%	excitations given T1,T2,TR,TE in ms.  dfreq is the resonant
%	frequency in Hz.  flip is in radians.
%

function [Msig,Mss] = sssignal(flip,T1,T2,TE,TR,dfreq)

Rflip = yrot(flip);
[Atr,Btr] = freeprecess(TR-TE,T1,T2,dfreq);
[Ate,Bte] = freeprecess(TE,T1,T2,dfreq);

% Let 	M1 be the magnetization just before the tip.
%	M2 be just after the tip.
%	M3 be at TE.
%
% then
%	M2 = Rflip * M1
%	M3 = Ate * M2 + Bte
%	M1 = Atr * M3 + Btr
%
% Solve for M3...
%
%	M3 = Ate*Rflip*Atr*M3 + (Ate*Rflip*Btr+Bte)

Mss = inv(eye(3)-Ate*Rflip*Atr) * (Ate*Rflip*Btr+Bte);
Msig = Mss(1)+i*Mss(2);

end