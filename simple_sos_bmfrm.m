function [img] = simple_sos_bmfrm(rf_data,fs,sos_t,sos_m,rx_delay,fnum,pos_trans,pos_z,pos_x,ang_t)

% Simple delay-and-sum beamformer which accounts for difference between assumed
% probe speed-of-sound and the speed-of-sound of the imaged medium,
% associated with the paper "Real-Time Speed-of-Sound Estimation In Vivo
% via Steered Plane Wave Ultrasound" in IEEE T-UFFC
% ver. 1.0  (Mar 2024 -- by Di Xiao, Pat De la Torre & Alfred Yu)
%
% simple_sos_bmfrm.m
% 
% Usage: 
% This is a function that takes in preprocessed RF data alongside the
% beamforming parameters and outputs an image
% 
% Syntax:
%   img = simple_sos_bmfrm(rf_data,fs,sos_t,sos_m,rx_delay,fnum,pos_trans,pos_z,pos_x,ang_t)
%
% Input:
%   rf_data         -- 2D array of postprocessed RF data (samples x channels)
%   fs              -- Sampling frequency of RF data
%   sos_t           -- Speed-of-sound used for transmission of steered plane wave
%   sos_m           -- Assumed speed-of-sound of the imaged medium
%   rx_delay        -- Receive delay corresponding to time of first element firing 
%   fnum            -- F-number used for aperture selection during beamforming
%   pos_trans       -- Position of linear array elements (channels x 1 array)
%   pos_z           -- Axial positions of pixels to be beamformed (Nz x 1 array)
%   pos_x           -- Lateral positions of pixels to be beamformed (Nx x 1 array)
%   ang_t           -- Intended steering angle for plane wave transmission
%
% Output:
%   img             -- 2D intensity image
%

% Automatically identifying the lengths of the input arrays
Nz = length(pos_z);
Nx = length(pos_x);
Nc = size(rf_data,2);
Ns = size(rf_data,1);

img = zeros(Nz,Nx);

% Changing which element is first to activate for pos/neg PW steering
if ang_t<0  
    wave_source = pos_trans(end);
else
    wave_source = pos_trans(1);
end

% Adjusting the transmitted angle through the medium based on Snell's law (Eq. (1) in paper)
ang_m = asind(sos_m*sind(ang_t)/sos_t); 

% Looping through all the pixels to be beamformed
for z = 1:Nz
    % Setting aperture based on selected F-number
    a = pos_z(z)/(2*fnum);
    for x = 1:Nx
        % Checking for whether the pixel is insonified by the plane wave
        insonify_check = 0;
        if ang_m <= 0
            if pos_x(x) < wave_source+pos_z(z)*tand(ang_m)
                insonify_check = 1;
            end
        else
            if pos_x(x) > wave_source+pos_z(z)*tand(ang_m)
                insonify_check = 1;
            end
        end
        
        % Beamforming for a single insonified pixel with sample linear interpolation
        if insonify_check == 1
            % Calculating the transmission time based on Eq. (3) in the associated paper
            tx_time = pos_z(z)*cosd(ang_m)/sos_m + (pos_x(x)-wave_source)*sind(ang_m)/sos_m;

            % Calculation of receive time based on geometric principles
            rx_time = sqrt(pos_z(z)^2 + (pos_x(x) - pos_trans).^2)/sos_m;

            % Calculation of best sample alongside linear interpolation factor
            best_samp = fs*(tx_time + rx_time - rx_delay);
            s_bot = floor(best_samp);
            s_interp = best_samp-s_bot;

            % Performing summation across channels for appropriately delayed and interpolated samples
            for c = 1:Nc
                if abs(pos_trans(c)-pos_x(x)) < a
                    img(z,x) = img(z,x) +  (1-s_interp(c))*rf_data(s_bot(c),c) + (s_interp(c))*rf_data(s_bot(c)+1,c);
                end
            end
        end
    end
end