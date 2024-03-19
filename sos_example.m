% Global Speed-of-Sound Estimation Example
% ver. 1.0  (Mar 2024 -- by Di Xiao, Pat De la Torre & Alfred Yu)
%
% sos_example.m
% 
% Usage: 
% The main script for estimating the speed-of-sound of roughly homogeneous
% media, associated with the paper "Real-Time Speed-of-Sound Estimation In Vivo
% via Steered Plane Wave Ultrasound" in IEEE T-UFFC.

%% Generating imaging parameters
gen_params

%% Processing SoS for each agar phantom
% Initializing some arrays and figure
load_fold = 'Data\';
loss_fig = figure;
sos_est_array = zeros(15,1);
sos_ref_array = zeros(15,1);

% Looping through the 15 agar phantoms
for idx = 1:15
    % Loading and preprocessing the data (Fig. 3, Step 1)
    disp(['Phantom ' num2str(idx)])
    load([load_fold 'Phantom' num2str(idx) '.mat'])
    filt_rf = preprocess_rf(rf_data,filt_coeff);

    % Beamforming the data for all transmit angles for all SoS candidates
    sos_img_tensor = zeros(Nz,Nx,length(tx_angles),length(sos_coarse));
    for sos_idx = 1:length(sos_coarse)
        disp(['Beamforming @ ' num2str(sos_coarse(sos_idx)) ' m/s'])

        % Adjusting the axial pixel positions for different SoS values (based on Eq. (4) in the paper)
        pos_z = pos_z_1540*sos_coarse(sos_idx)/sos_probe;

        % Performing the beamforming (Fig. 3, Step 2)
        for ang_idx = 1:length(tx_angles)
            sos_img_tensor(:,:,ang_idx,sos_idx) = simple_sos_bmfrm(filt_rf(:,:,ang_idx),fs,sos_probe,...
                sos_coarse(sos_idx),rx_delay,fnum,...
                pos_trans,pos_z,pos_x,tx_angles(ang_idx));
        end
    end

    % Calculating the inter-steering loss based on all the beamformed images (Fig. 3, Step 3)
    inter_steer_loss = calc_intersteering_loss(sos_img_tensor);

    % Fitting smoothing spline to beamformed SoS candidates (Fig. 3, Step 4)
    fitted_loss = fit(sos_coarse',inter_steer_loss,'smoothingspline','Normalize','on','SmoothingParam',1);

    % Estimating SoS based on minimum of fitted function (Fig. 3, Step 4)
    [~,I] = min(fitted_loss(sos_fine));
    sos_est = sos_fine(I);
    disp(['SoS Estimate: ' num2str(sos_est) ' m/s      SoS Reference: ' num2str(sos_meas) ' m/s'])

    % Plotting and storing values
    plot_loss(loss_fig,sos_coarse,sos_fine,inter_steer_loss,fitted_loss(sos_fine),sos_est,sos_meas)
    sos_est_array(idx) = sos_est;
    sos_ref_array(idx) = sos_meas;
end


plot_comparison(sos_est_array,sos_ref_array)
mean(sos_est_array-sos_ref_array)

%% Demonstration of refraction-modelling beamformer with CIRS phantom data

load('Data\CIRS.mat')

filt_rf = preprocess_rf(rf_data,filt_coeff);

% Beamforming the data for all tx angles with modified medium SoS
sos_img_tensor = zeros(Nz,Nx,length(tx_angles),length(sos_coarse));
for sos_idx = 1:length(sos_coarse)
    disp(['Beamforming @ ' num2str(sos_coarse(sos_idx)) ' m/s'])
    pos_z = pos_z_1540*sos_coarse(sos_idx)/sos_probe;
    for ang_idx = 1:length(tx_angles)
        sos_img_tensor(:,:,ang_idx,sos_idx) = simple_sos_bmfrm(filt_rf(:,:,ang_idx),fs,sos_probe,...
                                                                    sos_coarse(sos_idx),rx_delay,fnum,...
                                                                    pos_trans,pos_z,pos_x,tx_angles(ang_idx));
    end
end

figure
for sos_idx = 1:length(sos_coarse)
    for ang_idx = 1:length(tx_angles)
        img = sos_img_tensor(:,:,ang_idx,sos_idx);
        img = 20*log10(abs(img));
        imagesc(img, [max(img(:))-40, max(img(:))])
        colormap gray
        title(num2str(sos_coarse(sos_idx)))
        drawnow
    end
end

%% Defining Local Functions

% Local funtion for preprocessing the raw RF data through bandpass filtering and Hilbert transform
function rf = preprocess_rf(raw_rf,filt_coeff)
rf = hilbert(filtfilt(filt_coeff,1,raw_rf));
end

% Local function to calculate the inter-steering loss
function raw_loss = calc_intersteering_loss(img_tensor)

% Masking the pixels that are insonified by all transmissions
mask = ones(size(img_tensor));
mask(img_tensor == 0) = 0;
mask = prod(mask,3);
img_tensor = img_tensor.*mask;
img_tensor(img_tensor == 0) = NaN;

% Setting a normalization factor for each SoS candidate based on the number
% of pixels that are actually insonified for each adjusted angle
norm_factor = squeeze(sum(mask,[1 2]));

% Calculating intersteering loss based only on the pixels that are fully
% insonified for all tranmissions
coeff_var = std(img_tensor,0,3,'omitnan')./mean(abs(img_tensor),3,'omitnan');
raw_loss = squeeze(sum(coeff_var,[1 2],'omitnan'))./norm_factor;
end

% Local function for plotting the loss function for a given phantom
function plot_loss(fig_handle,sos_coarse,sos_fine,raw_loss,fit_loss,sos_est,sos_ref)
figure(fig_handle)
clf(fig_handle)
plot(sos_fine,fit_loss,'Linewidth',2)
hold on
scatter(sos_coarse,raw_loss,'Linewidth',2)
xline(sos_ref,':','Linewidth',2)
title(['Estimated: ' num2str(sos_est), ' m/s  Reference: ' num2str(sos_ref) , ' m/s'])
legend('Smoothing','Loss','Reference')
xlim([min(sos_coarse)-10 max(sos_coarse)+10]);
drawnow
end

% Local function for plotting the overall accuracy of SoS estimation method
function plot_comparison(sos_est_array,sos_ref_array)
figure
grid on
box on

scatter(sos_ref_array,sos_est_array,50,[0.6350 0.0780 0.1840],'^','filled','LineWidth',2);
axis equal tight

xlim([1500 1700]);
ylim([1500 1700]);
xticks(1500:25:1700)
yticks(1500:25:1700)
xlabel('SoS Reference Measurement')
ylabel('SoS Estimation')

hline = refline(1,0);
hline.LineStyle = '--';
hline.LineWidth = 2;
hline.Color = [0 0 0 0.25];
alpha(hline,.5)
grid on
drawnow
end