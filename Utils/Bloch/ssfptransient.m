%	== Example of transient SSFP response calculation using bloch.m ==
%	
clear all;close all;clc
%	--- Setup parameters ---
TR = .003;		% Sec.
Trf = 0.0005;		% 100 us "hard" RF pulse.
alpha = 30;		% Degrees.
gamma = 4258;		% Hz/G.
T1 = 1;			% Sec.
T2 = .1;		% Sec.
freq = [-200:1:200];	% Hz
N = 1000;
Tpad = (TR-Trf)/2;	% Sec.

% 	--- Setup B1, just delta-function RF ---
%
%	-- Time intervals are [Tpad, Trf, Tpad]
t = [Tpad Trf Tpad];

%	-- b1 is non-zero during Trf.
b1 = [0 pi/180*alpha/Trf/gamma/2/pi 0];	% Gauss.
b1_c = b1;

% === Calculate steady state, for comparison.
[mxss,myss,mzss] = bloch1(b1,0*b1,t,T1,T2,freq,0,1);
mss = mxss + 1j*myss;

% === Start with alpha/2, from equilibrium. ===
[mx,my,mz] = bloch1(+i*max(b1)/2,0,Trf,T1,T2,freq,0,0);
[mx_c,my_c,mz_c] = bloch1(+i*max(b1)/2,0,Trf,T1,T2,freq,mx,my,mz);

% Add a gradient induced zeroth order phase correction model
ga = 700000;
tx = 0.06;
ty = 0.06;
rad_ang = (pi / (((1 + sqrt(5)) / 2) + ga - 1));
model_girf = @(phi,theta_x,theta_y)(theta_x * cos(phi) + theta_y * sin(phi) );

for n=2:N
    b1 = b1 * exp(-1j * pi);
    b1_c = b1_c * exp(-1j * (pi + model_girf(rad_ang * n, tx,ty)));
	[mx(n,:),my(n,:),mz(n,:)] = bloch1(b1,0*b1,t,T1,T2,freq,0,0,mx(n - 1,:),my(n - 1,:),mz(n - 1,:));    
  	[mx_c(n,:),my_c(n,:),mz_c(n,:)] = bloch1(b1_c,0*b1_c,t,T1,T2,freq,0,0,mx_c(n - 1,:),my_c(n - 1,:),mz_c(n - 1,:));    
end

sig = mx + 1j * my;
sig_c = mx_c + 1j * my_c;

figure,plot(abs(sig(end,:)),'k');hold on;
plot(real(sig(end,:)),'r');
plot(imag(sig(end,:)),'b');
plot(abs(sig_c(end,:)),'k--');
plot(real(sig_c(end,:)),'r--');
plot(imag(sig_c(end,:)),'b--');


