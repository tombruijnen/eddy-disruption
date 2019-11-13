function sig = bSSFP_sim(sim,ec_model)
%Simulate steady-state bSSFP sequence for one tissue type and a range off
% off-resonances.

% Steady state initialization (fixed gradients) 
b1_0 = [0 pi/180*sim.alpha / sim.Trf / sim.gamma /2/pi 0];
mx = 0; my = 0; mz = 1;
for n = 1 : max([sim.dummies 1])
    b1 = b1_0 * exp(1j * mod(n,2) * sim.rfphase);
    [mx,my,mz] = bloch1(b1,0*b1,sim.t,sim.T1,sim.T2,sim.freq,0,0,mx,my,mz);
    
    % apply extra dephasing - because of eddy current
    tmp = (mx + 1j * my) * exp(1j * ec_model(sim.pe(1,:),sim.phi));
    mx = real(tmp);
    my = imag(tmp);
    
    if n > 1
        % Apply the dephasing from previous TR 
        tmp = (mx + 1j * my) * exp(-1j * ec_model(sim.pe(1,:),sim.phi));
        mx = real(tmp);
        my = imag(tmp);
    end
end

% Iterate over RF pulses
sig(1,:) = mx + 1j * my;
for n = 2 : sim.N 
    % apply RF-pulse
    b1_0 = [0 pi/180*sim.alpha / sim.Trf / sim.gamma /2/pi 0];	
    b1 = b1_0 * exp(1j * mod(n,2) * sim.rfphase);    
    [mx(2,:),my(2,:),mz(2,:)] = bloch1(b1,0*b1,sim.t,sim.T1,sim.T2,sim.freq,0,0,mx(1,:),my(1,:),mz(1,:));    
    
    % apply extra dephasing - because of eddy current    
    tmp = (mx(2,:) + 1j * my(2,:)) * exp(1j * ec_model(sim.pe(n,:),sim.phi));
    mx(2,:) = real(tmp);
    my(2,:) = imag(tmp);
    
    % Apply the dephasing from previous TR 
    tmp = (mx(2,:) + 1j * my(2,:)) * exp(-1j * ec_model(sim.pe(n - 1,:),sim.phi));
    mx(2,:) = real(tmp);
    my(2,:) = imag(tmp);
    
    sig(n,:) = mx(2,:) + 1j * my(2,:);
    mx(1,:) = [];
    my(1,:) = [];
    mz(1,:) = [];
end

% Demodulate the phase alternation
if sim.rfphase == pi
    for n = 1 : size(sig,1)
        if mod(n,2) > 0
            sig(n,:) = -sig(n,:);
        end
    end
end

% END
end

