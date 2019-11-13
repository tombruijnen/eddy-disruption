% Adjustable
sim.gamma   = 4258;                         % [Hz/G]
sim.TR      = G_t(end);                     % repetition time [s]
sim.Trf     = 0.0005;                       % rf pulse duration [s] 
sim.T1      = 1;                            % T1 value of spin [s]
sim.T2      = .1;                           % T2 value of spin [s]
sim.dfreq   = 1;                            % Spectral resolution [Hz]
sim.freq    = -400:sim.dfreq:400;           % Spectral range [Hz]
%sim.freq    = 0;
sim.N       = 512;                          % Number of rf pulses [-]
sim.alpha   = 35;                           % Flip angle [degrees]
sim.dummies = 0;                          % Dummy pulses [-]
sim.rfphase = pi;                           % Phase cycle scheme [rad]
sim.idim    = 128;                          % Image dimension [-]

% Calculate from above
sim.Tpad    = (sim.TR - sim.Trf)/2;        
sim.t       = [sim.Tpad sim.Trf sim.Tpad];
sim.pe      = ones(sim.N,1);
