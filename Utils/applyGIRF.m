function [cwf,b0_ec,ph_ec] = applyGIRF(t,nom,girf)
%Apply girfs on nominal waveform and return corrected waveform and B0
% modulations.
% Everything is in seconds and T/m/s
%    
% Version: 20171110 
% Author: Tom Bruijnen
% Contact: t.bruijnen@umcutrecht.nl

% check input

% Zeropad input with 30 ms on both sides to a 1 us grid
t_pad = 50E-03;
t_zp = (0:1E-06:2 * t_pad + t(end))';
nom_zp = interp1(t + t_pad,nom,t_zp,'phchip','extrap');
dt = abs(t_zp(2) - t_zp(1));

% Fourier transform the zeropadded waveform
F_nom=ifftshift(fft(fftshift(nom_zp,1),numel(t_zp),1),1) / numel(t_zp);

% Generate frequency vector of nominal waveform
df = (1/(dt)) / numel(t_zp);
f_nom = df*(0:numel(t_zp)-1);
f_nom = f_nom - df*ceil((numel(t_zp)-1)/2); 

% Resample the GIRFs to the nominal waveform frequencies
for ax=1:3
    g0(:,ax)=interp1(girf.freq,girf.girf0(:,ax),f_nom);g0(isnan(g0))=0;
    g1(:,ax)=interp1(girf.freq,girf.girf1(:,ax),f_nom);g1(isnan(g1))=0;
end

% Apply GIRFs
for ax1 = 1 : 3 % axis
    for ax2 = 1 : 3 % girfs
        F_b0(:,ax2,ax1) = F_nom(:,ax1).*g0(:,ax2);
        F_cwf(:,ax2,ax1) = F_nom(:,ax1).*g1(:,ax2);
    end
end

% Translate zeroth order to B0 modulations and to phase errors
b0_ec_zp = real(fftshift(ifft(fftshift(F_b0,1),numel(t_zp),1),1)) * numel(t_zp);
ph_ec_zp = 2 * pi * dt * cumsum(b0_ec_zp,1);

% Apply first order to gradient waveform
cwf_zp=real(fftshift(ifft(fftshift(F_cwf,1),numel(t_zp),1),1)) * numel(t_zp);

% Remove zero paddings
for ax1 = 1 : 3
    for ax2 = 1 : 3
        cwf(:,ax1,ax2) = interp1(t_zp - t_pad,cwf_zp(:,ax1,ax2),t);
        b0_ec(:,ax1,ax2) = interp1(t_zp - t_pad,b0_ec_zp(:,ax1,ax2),t);
        ph_ec(:,ax1,ax2) = interp1(t_zp - t_pad,ph_ec_zp(:,ax1,ax2),t);
    end
end

% END
end