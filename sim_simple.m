%% ------------------
%% FMCW RX Simulation
%% ------------------
fminR = 17e3;
B = 5e3;
Fs = 48000;
vs = 340;
sampleInterval=0.030; % 30 ms
simLenS =  10;

radius = 0.05;  % array radius (m)
Nr = 8; % 8 microphones
rtDist = 8; % approx. round-trip distance from sound source to microphones
xPosWall = 1; % side wall x position for multipath
mpAngle = 60; % multipath angle in degree
addNoise = true;

fmaxR = fminR + B;
fc = (fminR + fmaxR)/2;
Ts=1/Fs;
K=sampleInterval/Ts;

t=(0:K-1)*Ts;
%nChirps = floor(simLenS * Fs / K);
nChirps = 5; %1; %2;
simLenS = nChirps * K / Fs;

St0 = cos(2*pi*(1/2*t.^2*B/sampleInterval+fminR*t)); % original Tx signal
St = repmat(St0, 1 , nChirps);


[m_xPos, m_yPos, m_zPos, rxarray]  = generate_rx_tx_positions(Nr, radius, fmaxR); % for now generates circular array
[Sr_noise, Sr, s_Pos, distance] = generate_rx_fmcw(fminR, B, Fs, vs, sampleInterval, nChirps, Nr, rxarray, rtDist, xPosWall, mpAngle, true, false);

%[Sr_noise1, Sr1, s_Pos1, distance1] = generate_rx_fmcw1(fminR, B, Fs, vs, sampleInterval, nChirps, Nr, rxarray, rtDist, xPosWall, mpAngle);
%[xx, ww] = opt_beam2(Sr1);

% figure;
% plot(Sr_noise(1, :));
% figure; spectrogram(Sr_noise(1,:),'yaxis',128,120,128,Fs)

% Apply fft filter
for mic = 1 : Nr
    Sr_noise(mic, :) = fftFilter(Sr_noise(mic, :),Fs,fminR,fmaxR,50);
end

%figure; spectrogram(Sr_noise(1,:),'yaxis',128,120,128,Fs)

% % de-chirping
[f, profile_mic1] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, Sr_noise(1,:), 'No Beam at Mic1');
%[f, profile_mic2] = dechirp_fmcw(Fs, fminR, B, vs, sampleInterval, nChirps, Sr_noise(2,:));


%incidentAz = 90;
incidentAz = 90; 

[y_DAS, y_MVDR, y_LCMV, y_LP, y_FR] = beamform(incidentAz, fc, vs, Fs, Sr_noise, rxarray, fmaxR, m_xPos, m_yPos, m_zPos, Nr);

w = opt_beam(Sr_noise); % my own optimization code just takes in the collected signal
y_MINE = w * Sr_noise;

[f, profile_DAS] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_DAS, 'DAS');
[f, profile_MVDR] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_MVDR, 'MVDR');
[f, profile_LCMV] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_LCMV, 'LCMV');
[f, profile_LP] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_LP, 'LP');
[f, profile_MINE] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_MINE, 'MINE');
%[f, profile_FR] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_FR, 'FR');
%[f, profile_PS] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_PS, 'PS');

figure;
dist = vs*f*sampleInterval*1000/B;
plot(dist, profile_mic1(1,:));
hold on;
plot(dist, profile_DAS(1,:));
plot(dist, profile_MVDR(1,:));
plot(dist, profile_LCMV(1,:));
plot(dist, profile_LP(1,:));
plot(dist, profile_MINE(1,:));
%plot(dist, profile_FR(1,:));
%plot(dist, profile_PS(1,:));
xlim([0 rtDist*2]);
title ('FMCW Profile')
xlabel('Distance (m)')
ylabel('Amplitude')
legend("No Beam at Mic1","DAS","MVDR", "LCMV", "LP", "MINE") %"FR") %, "PS")

% %% compare beamforming performance varing num mics used from 2 - Nr
% profile_DAS_mics = [];
% % profile_MVDR_mics = [];
% % profile_LCMV_mics = [];
% profile_LP_mics = [];
% %profile_PS_mics = [];
%
% for num_mics = 2:Nr
%     [y_DAS1, y_MVDR1, y_LCMV1, y_LP1, y_FR1] = beamform(incidentAz, fc, vs, Fs, Sr_noise, rxarray, fmaxR , m_xPos, m_yPos, m_zPos, num_mics);
%     
%     [f, profile_DAS1] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_DAS1, strcat('DAS w ', num2str(num_mics), ' mics'));
% %     [f, profile_MVDR1] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_MVDR1, 'MVDR');
% %     [f, profile_LCMV1] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_LCMV1, 'LCMV');
%      [f, profile_LP1] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_LP1, strcat('LP w ', num2str(num_mics), ' mics'));
% %     [f, profile_PS1] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_PS1, 'PS');
%  
%     profile_DAS_mics = [profile_DAS_mics; profile_DAS1(1,:)];
%     profile_LP_mics = [profile_LP_mics; profile_LP1(1,:)];
% end
% 
% figure;
% dist = vs*f*sampleInterval*1000/B;
% plot(dist, profile_mic1(1,:));
% hold on;
% for num_mics = 2:Nr
%     if(mod(num_mics, 2) == 0) 
%         plot(dist, profile_DAS_mics(num_mics-1,:));
%     end
% end 
% title ('FMCW Profile w varying num mics beamforming (DAS)')
% xlim([0 rtDist*2]);
% xlabel('Distance (m)')
% ylabel('Amplitude')
% %legend("No beam", "2 Mics","3 Mics","4 Mics", "5 Mics", "6 Mics", "7 Mics", "8 Mics")
% legend("No beam", "2 Mics","4 Mics", "6 Mics", "8 Mics")
% 
% 
% figure;
% dist = vs*f*sampleInterval*1000/B;
% plot(dist, profile_mic1(1,:));
% hold on;
% for num_mics = 2:Nr
%     if(mod(num_mics, 2) == 0) 
%         plot(dist, profile_LP_mics(num_mics-1,:));
%     end
% end 
% title ('FMCW Profile w varying num mics beamforming (LP)')
% xlim([0 rtDist*2]);
% xlabel('Distance (m)')
% ylabel('Amplitude')
% %legend("No beam", "2 Mics","3 Mics","4 Mics", "5 Mics", "6 Mics", "7 Mics", "8 Mics")
% legend("No beam", "2 Mics","4 Mics", "6 Mics", "8 Mics")
