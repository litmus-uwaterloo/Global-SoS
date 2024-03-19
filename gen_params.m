% Parameter Generation for Global Speed-of-Sound Estimation Example
% ver. 1.0  (Mar 2024 -- by Di Xiao, Pat De la Torre & Alfred Yu)
%
% gen_params.m
% 
% Usage: 
% The helper script for setting the parameters for speed-of-sound estimation, 
% associated with the paper "Real-Time Speed-of-Sound Estimation In Vivo
% via Steered Plane Wave Ultrasound" in IEEE T-UFFC.

%%
clear

%% Setting initial parameters

% Transmission settings
tx_angles = [-15:3:-3 3:3:15]; % steering transmission angles in degrees
fc = 5e6; % transmit center frequency
pitch = 0.3048e-3; % probe element pitch
sos_probe = 1540; % assumed SoS of probe for transmission

% Receive parameters
fs = 40e6; % sampling frequency of ADC
rx_delay = 4.1e-6; % receive delay was characterized experimentally in water

% Data sizes
Nz = 128; % axial pixels
Nx = 128; % lateral pixels
Nc = 128; % number of transducer channels
Ns = 3120; % number of samples per channel

% Postprocessing and beamforming parameters
pos_trans = pitch*linspace(-(Nc-1)/2,(Nc-1)/2,Nc); % transducer lateral positions (linear array)
filt_coeff = fir1(40,[3.5e6 6.5e6]/(fs/2)); % filter coefficients for bandpass filtering
sos_coarse = 1400:20:1740; % SoS candidates for beamforming (coarse SoS grid)
sos_fine = 1400:0.1:1740; % SoS values used for finding the minimum of the fitted SoS loss
pos_z_1540 = linspace(20e-3, 45e-3, Nz); % Axial pixel positions used for beamforming (at 1540 m/s)
pos_x = linspace(-18e-3, 18e-3, Nx); % Lateral pixel positions used for beamforming
fnum = 1.4; % F-number to set transmit aperture