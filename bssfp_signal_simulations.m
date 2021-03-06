%% bSSFP signal evolutions using Bloch simulations and GIRFs
% SUMMARY: The gradient impulse response function (GIRF) can be used to
% predict additional phase accumulation during the pulse sequence, which can
% subsequently be included in Bloch simulations. This additional phase
% accumulation affects can considerable affect the convergence of spins to
% the steady-state. Some of these effects are studied in this script.
% Figures similar to figure 4 from the article are produced.
%
% NOTES:     
%             
%
% Acronyms: 
%   - LIN = linear phase encodes
%   - RND = random phase encodes
%   - GA  = golden angle
%   - RAD = radial
%  
%
% Contact: T.Bruijnen@umcutrecht.nl | University Medical Center Utrecht
% Department of Radiotherapy, University Medical Center Utrecht, Utrecht, the Netherlands
% Computational Imaging Group for MRI diagnostics and therapy, Centre for
% Image Sciences, University Medical Center Utrecht, Utrecht, The Netherlands
%
% Init: T.Bruijnen - 20191113

addpath(genpath(pwd))
load('Data/labels_simulation.mat')

% Get default simulation parameters - these are all adjustable
run get_simulation_parameters;

% Use GIRF to calculate phase errors over the TR
% Note that this shows the errors for the maximum gradient (i.e. max PE)
[~,~,d_phi] = applyGIRF(G_t,G,girf);

% Interpolate the phase error at the end of the TR - equation (3) of paper
% Note that we pick one gradient axis now to function as phase encode axis
d_phi_max = d_phi(end,1:2,2);

% Define eddy current model 
% Note that this is linear scaling of the phase error w.r.t phase encode
ec_model = @(scale,phi)(sum(scale .* phi));

% Define the sampling schemes
samp_lin = matrix_to_vec(linspace(-1,1,sim.N)');
samp_rnd = matrix_to_vec(demax(unique(ceil(sim.N*rand(10000,1)),'stable')-(sim.N/2 + 1)));
samp_rad = radial_trajectory([2*sim.N sim.N],1); % 1 = first golden angle
samp_rad = demax(squeeze(samp_rad(1:2,1,:)))';

% Run the simulation without eddy current
sim.phi = 0;
sim.pe  = samp_lin;
signal = bSSFP_sim(sim,ec_model);

% Run the simulation with eddy currents & LIN
sim.phi = d_phi_max(1);
signal_lin = bSSFP_sim(sim,ec_model);

% Run the simulation with eddy currents & RND
sim.pe  = samp_rnd;
signal_rnd = bSSFP_sim(sim,ec_model);

% Run the simulation with eddy currents & GA-RAD
sim.pe  = samp_rad;
signal_rad = bSSFP_sim(sim,ec_model);

%% Visualization
run visualization_bssfp_signal_simulation.m
