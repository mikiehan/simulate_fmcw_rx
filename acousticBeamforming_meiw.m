% parameters
c = 340;        % sound speed
fs = 48000;     % sampling rate
duration = 5;   % duration time
Nr = 8;         % # of antenna element of Rx
Ns = 1;         % # source signal
incidentAz = 90; % incident angle (from where sound source is coming)

useMatrixVoice = true; % if false UCA is used
useTrace = false; % if false generated FMCW rx signal is used
addNoise = false;
if useTrace % when using collected trace it is matrixVoice
    useMatrixVoice = true;
end

fmin = 10e3;
fmax = 15e3;
B = fmax - fmin;
chirp_dur = 0.030; % 0.030 ms
fc = (fmin + fmax)/2;

lambda = c/fc;  % wavelength

%addNoise = true;

%%
% Define a Uniform Microphone Array
microphone = phased.OmnidirectionalMicrophoneElement('FrequencyRange',[20 fmax]);

if(useMatrixVoice)
    
    m_xPos = [0, -38.13, -20.98, 11.97, 35.91, 32.81, 5.00, -26.57] ./ 1000;
    m_yPos = [0, 3.58, 32.04, 36.38, 13.32, -19.77, -37.97, -27.58] ./ 1000;
    m_zPos = zeros(1, Nr);
    [m_gamma, m_l] = cart2pol(m_xPos, m_yPos);
    nvec = [90, m_gamma(2:end)./pi * 180];
    
    rxarray = phased.ConformalArray(...
        'ElementPosition',[m_xPos;m_yPos;m_zPos],...
        'ElementNormal',[nvec;zeros(1,Nr)],...
        'Element', microphone);
    
else
    radius = 0.02;  % array radius (m)
    rxarray = phased.UCA('NumElements', Nr ,'Radius', radius, "Element", microphone);
    pos = getElementPosition(rxarray);
    m_xPos = pos(1,:); % x-coord
    m_yPos = pos(2,:); % y-coord
    m_zPos = pos(3,:); % z-coord
    
end

figure(1);
viewArray(rxarray, 'ShowNormals', true);


% Define object position
s_Pos = [0, 0.9144, 0]; % source position in meter

% plot setup
figure(2);
hold on;
for ni = 1:Nr
    scatter(m_xPos(ni),m_yPos(ni),'*','b')
    text(m_xPos(ni),m_yPos(ni),num2str(ni))
end
scatter(s_Pos(1),s_Pos(2),'d')
text(s_Pos(1),s_Pos(2),string(s_Pos(2)))
xlim([-0.5 0.5])
ylim([-0.5 1])
hold off
title('array setup')
xlabel('x(m)')
ylabel('y(m)')


t = 0:1/fs:duration;
t = t(1:end-1)';

if useTrace
    trial = 1;
    trimSec = 0;
    Fs = 48000;
    %inch = 72; % [48, 60, 72, 84, 96, 108, 120, 132, 144]; % for 18-20kHz
    inch = 36; %36; % [36, 48, 60, 72, 84, 96];
    prefix = sprintf('fmcw_16k_20k_30_%dinch_home_v_%d', inch, trial);
    prefix = 'fmcw_10k_15k_25_36inch_room_samedev_1_48000'; %_channel_7_cleaned'
    
    display(prefix);
    dir = strcat('/Users/profhan/Research/fmcw_aoa_codes/aoa_beamforming_script/08-09-new/home_wav/');
    dir = '~/Downloads/';
    display(strcat(dir , prefix , '_channel_0.wav'));
    [y0,Fs] = audioread(strcat(dir , prefix , '_channel_0_cleaned.wav'));
    [y1,Fs] = audioread(strcat(dir , prefix , '_channel_1_cleaned.wav'));
    [y2,Fs] = audioread(strcat(dir , prefix , '_channel_2_cleaned.wav'));
    [y3,Fs] = audioread(strcat(dir , prefix , '_channel_3_cleaned.wav'));
    [y4,Fs] = audioread(strcat(dir , prefix , '_channel_4_cleaned.wav'));
    [y5,Fs] = audioread(strcat(dir , prefix , '_channel_5_cleaned.wav'));
    [y6,Fs] = audioread(strcat(dir , prefix , '_channel_6_cleaned.wav'));
    [y7,Fs] = audioread(strcat(dir , prefix , '_channel_7_cleaned.wav'));
    
    y0 = y0(trimSec*Fs+1:(trimSec+duration)*Fs);
    y1 = y1(trimSec*Fs+1:(trimSec+duration)*Fs);
    y2 = y2(trimSec*Fs+1:(trimSec+duration)*Fs);
    y3 = y3(trimSec*Fs+1:(trimSec+duration)*Fs);
    y4 = y4(trimSec*Fs+1:(trimSec+duration)*Fs);
    y5 = y5(trimSec*Fs+1:(trimSec+duration)*Fs);
    y6 = y6(trimSec*Fs+1:(trimSec+duration)*Fs);
    y7 = y7(trimSec*Fs+1:(trimSec+duration)*Fs);
    
    rx = [y0 y1 y2 y3 y4 y5 y6 y7];
    noise = ones(size(rx))*1e-10;
else
    dir = '~/Downloads/';
    prefix =strcat('fmcw_sim_fmin', num2str(fmin), '_B' , num2str(B), '_', num2str(chirp_dur), '_', num2str(fs));
    rx_raw = gen_FMCW_rx_sim(s_Pos, m_xPos, m_yPos, fmin, B, chirp_dur, fs, duration);
    t = t(1:length(rx_raw));
    
    if(addNoise)
        noise = 0.5 * randn(size(rx_raw));
    else
        noise = ones(size(rx_raw))*1e-10;
    end
    
    rx = rx_raw + noise;
end

%%
% convert to complex data
x = zeros(size(rx));
%x_raw = zeros(size(rx_raw));
for mic = 1 : Nr
    x(:,mic) = fftFilterSingleSide(rx(:,mic)',fs,fmin,fmax,100);
    %    x_raw(:,mic) = fftFilterSingleSide(rx_raw(:,mic)',fs,fc-200,fc+200,100);
end
% x = rx;
%%
thetaScanningAngles = 90; % 0:1:360;
phiScanningAngles = 0:1:180; % 90
escan = steeringVector(m_xPos, m_yPos, m_zPos, fc, c, thetaScanningAngles, phiScanningAngles);

%%
incidentAngle = [incidentAz; 90]
% plot beamform patterns
ebi = steeringVector(m_xPos, m_yPos, m_zPos, fc, c, incidentAngle(1), incidentAngle(2));
ebi = reshape(ebi, 1, 1, 8);


beamformer = phased.PhaseShiftBeamformer('SensorArray',rxarray,...
    'OperatingFrequency',fc,'PropagationSpeed',c,...
    'Direction',[incidentAz; 90],'WeightsOutputPort',true);
[y_PS,w_PS] = beamformer(x);


w_DAS = weightingVectorDAS(ebi);
w_MVDR = weightingVectorMVDR(x.', ebi);
w_LCMV = weightingVectorLCMV(x.');
w_LP = weightingVectorLP(x.');
AF_DAS = arrayFactor(w_DAS, escan);
AF_MVDR = arrayFactor(w_MVDR, escan);
AF_LCMV = arrayFactor(w_LCMV, escan);
AF_LP = arrayFactor(w_LP, escan);
AF_PS = arrayFactor(w_PS, escan);

figure(3);
% polarplot(deg2rad(thetaScanningAngles),AF_DAS')
% hold on
% polarplot(deg2rad(thetaScanningAngles),AF_MVDR')
% polarplot(deg2rad(thetaScanningAngles),AF_LCMV')
% polarplot(deg2rad(thetaScanningAngles),AF_LP')
% polarplot(deg2rad(thetaScanningAngles),AF_PS')
% hold off

polarplot(deg2rad(phiScanningAngles),AF_DAS')
hold on
polarplot(deg2rad(phiScanningAngles),AF_MVDR')
polarplot(deg2rad(phiScanningAngles),AF_LCMV')
polarplot(deg2rad(phiScanningAngles),AF_LP')
polarplot(deg2rad(phiScanningAngles),AF_PS')
hold off

%legend(["DAS", "MVDR"]);
legend(["DAS", "MVDR", "LCMV", "LP", "PS"])
%%
% Beamforming
y_DAS = x*squeeze(w_DAS);
y_MVDR = x*squeeze(w_MVDR);
y_LCMV = x*squeeze(w_LCMV);
y_LP = x*squeeze(w_LP);

% figure(4);
% plot(t,real(x(:,1)), t,real(x_raw(:,1)));
% xlim([0 0.001]);
% legend('Raw','Raw no noise')


figure(4);
plot(t,real(x(:,2)),'b:',t,real(y_DAS),'g', t, real(y_MVDR),'c', t, real(y_LCMV), 'y', t, real(y_LP), 'r');
xlim([0 0.005]);
xlabel('Time (s)')
ylabel('Amplitude')
legend('Raw','DAS','MVDR', "LCMV", "LP")

%%
% calculate array gain
snr_raw = snr(real(x(:,2)), noise(:,2))
snr_DAS = snr(real(y_DAS), noise(:,2))
snr_MVDR = snr(real(y_MVDR), noise(:,2))
snr_LCMV = snr(real(y_LCMV), noise(:,2))
snr_LP = snr(real(y_LP), noise(:,2))
snr_PS = snr(real(y_PS), noise(:,2))

disp(snr_DAS - snr_raw);
disp(snr_LCMV - snr_raw);
disp(snr_LP - snr_raw);

figure(5);
plot(t,real(x(:,2)),'b:',t,real(y_DAS), 'g', t, real(y_LCMV), 'c', t, real(y_LP), 'r');
xlim([0 0.001]);
xlabel('Time (s)')
ylabel('Amplitude')
legend('Raw','DAS', "LCMV", "LP")

% figure(5);
% plot(t,real(y_DAS));
% xlim([0 0.1]);
%
% figure(6);
% plot(t,real(x(:,1)));
% xlim([0 0.1]);
normalized_yDAS = real(y_DAS)./max(abs(real(y_DAS)));
normalized_yLCMV = real(y_LCMV)./max(abs(real(y_LCMV)));
normalized_yLP = real(y_LP)./max(abs(real(y_LP)));

figure; spectrogram(y_DAS,'yaxis',128,120,128,fs);
%audiowrite(strcat(dir , prefix , '_channel_8_das.wav'), real(y_DAS), fs);
audiowrite(strcat(dir , prefix , '_channel_8_das.wav'), normalized_yDAS, fs);

%audiowrite(strcat(dir , prefix , '_channel_8_lcmv.wav'), real(y_LCMV), fs);
audiowrite(strcat(dir , prefix , '_channel_8_lcmv.wav'), normalized_yLCMV, fs);
%audiowrite(strcat(dir , prefix , '_channel_8_lp.wav'), real(y_LP), fs);
audiowrite(strcat(dir , prefix , '_channel_8_lp.wav'), normalized_yLP, fs);

audiowrite(strcat(dir , prefix , '_channel_0.wav'), rx_raw(:,1), fs);
