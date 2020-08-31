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
sOrigin = 8; % approx. 5 meter away (round-trip) from sound source to microphones
addNoise = true;

fmaxR = fminR + B;
Ts=1/Fs;
K=sampleInterval/Ts;

t=(0:K-1)*Ts;
%nChirps = floor(simLenS * Fs / K);
nChirps = 2;
simLenS = nChirps * K / Fs;

St0 = cos(2*pi*(1/2*t.^2*B/sampleInterval+fminR*t)); % original Tx signal
St = repmat(St0, 1 , nChirps);


[s_Pos, m_xPos, m_yPos, m_zPos, rxarray, distance]  = generate_rx_tx_positions(Nr, radius, sOrigin, fmaxR); % for now generates circular array
[Sr_noise, Sr] = generate_rx_fmcw(fminR, B, Fs, vs, sampleInterval, nChirps, Nr, distance, false);

% figure;
% plot(Sr_noise(1, :));
% figure; spectrogram(Sr_noise(1,:),'yaxis',128,120,128,Fs)

% Apply fft filter 
for mic = 1 : Nr
    Sr_noise(mic, :) = fftFilter(Sr_noise(mic, :),Fs,fminR,fmaxR,50);
end

%figure; spectrogram(Sr_noise(1,:),'yaxis',128,120,128,Fs)

% % de-chirping
[f, profile_mic1] = dechirp_fmcw(sOrigin, Fs, fminR, B, vs, sampleInterval, nChirps, Sr_noise(1,:), 'No Beam at Mic1');
%[f, profile_mic2] = dechirp_fmcw(Fs, fminR, B, vs, sampleInterval, nChirps, Sr_noise(2,:));


incidentAz = 90;
fc = (fminR + fmaxR)/2;
[y_DAS, y_MVDR, y_LCMV, y_LP, y_PS] = beamform(incidentAz, fc, vs, Sr_noise, rxarray, m_xPos, m_yPos, m_zPos);

[f, profile_DAS] = dechirp_fmcw(sOrigin, Fs, fminR, B, vs, sampleInterval, nChirps, y_DAS, 'DAS');
[f, profile_MVDR] = dechirp_fmcw(sOrigin, Fs, fminR, B, vs, sampleInterval, nChirps, y_MVDR, 'MVDR');
[f, profile_LCMV] = dechirp_fmcw(sOrigin, Fs, fminR, B, vs, sampleInterval, nChirps, y_LCMV, 'LCMV');
[f, profile_LP] = dechirp_fmcw(sOrigin, Fs, fminR, B, vs, sampleInterval, nChirps, y_LP, 'LP');
[f, profile_PS] = dechirp_fmcw(sOrigin, Fs, fminR, B, vs, sampleInterval, nChirps, y_PS, 'PS');

figure;
dist = vs*f*sampleInterval*1000/B;
plot(dist, profile_mic1(1,:)); 
hold on;
plot(dist, profile_DAS(1,:));
plot(dist, profile_MVDR(1,:));
plot(dist, profile_LCMV(1,:));
plot(dist, profile_LP(1,:));
plot(dist, profile_PS(1,:));
xlim([0 sOrigin*2]);
title ('FMCW Profile')
xlabel('Distance (m)')
ylabel('Amplitude')
legend("No Beam at Mic1","DAS","MVDR", "LCMV", "LP", "PS")

