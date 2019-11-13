%% Process gradient impulse response measurement data
% SUMMARY: Process thin slice method measurements data from a 1.5T Philips
% scanner and reconstruct the zeroth and first order gradient impulse
% response functions (GIRFs).
%
% NOTES:    - Data has been coil compressed and averaged already. 
%             This was done to reduce the size of the data for upload.
%
%
% Contact: T.Bruijnen@umcutrecht.nl | University Medical Center Utrecht
% Department of Radiotherapy, University Medical Center Utrecht, Utrecht, the Netherlands
% Computational Imaging Group for MRI diagnostics and therapy, Centre for
% Image Sciences, University Medical Center Utrecht, Utrecht, The Netherlands
%
% Init: T.Bruijnen - 20190912 

%% Raw data processing to obtain zeroth and first order responses
% Load data - dimensions: [nsamples ngradients 6 axis]
load('Data/kdata_girf_measurement.mat')
load('Data/labels_girf_measurement.mat')

% Phase unwrapping
kdata = unwrap(angle(kdata));

% Duyn's method as described in:
% Jeff H. Duyn,Yihong Yang, Joseph A. Frank, and Jan Willem van der Veen
% Simple Correction Method for k-Space Trajectory Deviations in MRI

% Subtract background phase, p = positve and n = negative
slice_1_p = squeeze(kdata(:,:,1,:) - kdata(:,:,6,:));
slice_1_n = squeeze(kdata(:,:,2,:) - kdata(:,:,6,:));
slice_2_p = squeeze(kdata(:,:,3,:) - kdata(:,:,5,:));
slice_2_n = squeeze(kdata(:,:,4,:) - kdata(:,:,5,:));

% Average positive and negative waveforms
slice_1_duyn = (slice_1_p + -1*slice_1_n) / 2;
slice_2_duyn = (slice_2_p + -1*slice_2_n) / 2;
first_order  = (slice_1_duyn + -1*slice_2_duyn) / 2;

% Derivative and convert to phase in radians
first_order = cat(1,zeros(1,labels.ngradients,3), diff(first_order,1,1)/labels.OffCenter)...
    / ((labels.t_adc(2)-labels.t_adc(1))*(267.513*10^6));

% Brodsky's method as described in:
% Ethan K. Brodsky, Jessica L. Klaers, Alexey A. Samsonov, 
%   Richard Kijowski and Walter F. Block
% Rapid Measurement and Correction of Phase Errors from B0 Eddy Currents: 
%   Impact on Image Quality for Non-Cartesian Imaging

% Subtract positive and negative waveforms as described by:
slice_1 = slice_1_p - slice_1_n;
slice_2 = slice_2_p - slice_2_n;

% Subtract both slices and divide by 4
zeroth_order = squeeze((slice_1 + slice_2) / 4);

% Derivative and convert to phase in radians
zeroth_order = cat(1,zeros(1,labels.ngradients,3),diff(zeroth_order)...
    *1/(labels.t_adc(2)-labels.t_adc(1)) * 1 / (2 * pi)); % [Hz]

% Visualization
run visualization_kdata_processing.m

%% Gradient impulse response function fit
% Interpolate measured waveforms to the high resolution 
for ax = 1 : 3
    for g = 1 : labels.ngradients
        zeroth_order_us(:,g,ax) = interp1(labels.t_adc,zeroth_order(:,g,ax),labels.t,'linear','extrap');
        first_order_us(:,g,ax)  = interp1(labels.t_adc,first_order(:,g,ax),labels.t,'linear','extrap');        
    end
end
    
% Get frequency vector of the nominal waveform G
f         = (1/labels.t(end)) * (0:numel(labels.t)-1);
f         = f - f(numel(labels.t) / 2 + 1);
f_storage = linspace(-3E04,3E04,10001); f_storage(end) = []; % Frequency vector to store GIRFs

% Fourier transform measurements and nominal waveform G
freq_input   = ifftshift(fft(fftshift(labels.G,1),[],1),1) / sqrt(numel(labels.t));
freq_zeroth  = double(ifftshift(fft(fftshift(zeroth_order_us,1),[],1),1) / sqrt(numel(labels.t)));
freq_first   = double(ifftshift(fft(fftshift(first_order_us,1),[],1),1) / sqrt(numel(labels.t)));

% Fit girfs as described by Vannesjo et al.
% Signe J. Vannesjo,Maximilan Haeberlin,Lars Kasper,Matteo Pavan,
%   Bertram J. Wilm,Christoph Barmet,Klaas P. Pruessmann
% Gradient system characterization by impulse response measurements with a dynamic field camera
for ax = 1 : 3    
    % Equation (12) of the paper
    I          = freq_input;
    O_0        = freq_zeroth(:,:,ax);
    O_1        = freq_first(:,:,ax);
    tmp_zeroth = sum(conj(I) .* O_0,2) ./ sum(abs(I) .* abs(I),2);
    tmp_first  = sum(conj(I) .* O_1,2) ./ sum(abs(I) .* abs(I),2);
    
    % Store only relevant frequencies (e.g. -30kHz till 30kHz)
    girf_zeroth(:,ax) = interp1(f,tmp_zeroth,f_storage,'linear');
    girf_first(:,ax)  = interp1(f,tmp_first,f_storage,'linear');
end

% Visualization
run visualization_girf_calculation.m












