%% OVERVIEW :
%% This script generates direct path and multipath signal for FMCW %%
%% (only one multipath per mic)                                    %%
%% and performs receiver-side beamforming on                       %%
%% both the combined signal and the decomposed signal              %%  
%% (assumes perfect knoweldge on the decomposed signal)            %%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 0: Settings for FMCW Signal %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fminR = 17e3;
B = 5e3;
Fs = 48000;
vs = 340;
sampleInterval=0.030; % 30 ms
%simLenS =  10;

nChirps = 2;

Nr = 4; % 4 microphones
rtDist = 8; % approx. round-trip distance from sound source to microphones
wallDist = 1; % side wall x position for multipath
addNoise = false; % for now no noise
SNR_db = 5; % not used

fmaxR = fminR + B;
fc = (fminR + fmaxR)/2;

lambda = vs/fc;
spacing = lambda/2;

Ts=1/Fs;
K=sampleInterval/Ts;
simLenS = nChirps * K / Fs;

t=(0:K-1)*Ts;

% mic positions
m_xPos = zeros(1, Nr);
m_yPos = spacing * (0:Nr-1);
m_zPos = zeros(1, Nr);

% object of reflection (assume point object)
o_xPos = rtDist/2;
o_yPos = 0;
o_zPos = 0;

distance = rtDist/2;
%att = 1; % for now just 1


%% STEP 1: Generate direct path signal (x_direct) at 90 deg                %%
%% x_direct contains direct path signal at each mic for one chirp duration %%
%% Since the direct path has the same distance for each mic at 90 deg      %%
%% direct signal (x_direct) is the same at each mic                        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


att = distance^-4; % attenuation
%att = 1; %ignore attenuation
tau =  distance*2/vs;
t_tau = t-tau; % time stamp for each sample for direct path

x_direct = zeros(Nr, length(t)); % each row contains direct path signal at mic_row
for i=1:Nr
    %x(i, :) = att * sin(2*pi*fc*t_tau);
    %x(i, :) = att * exp(1i*2*pi*fc*t_tau);
    x_direct(i, :) = att * exp(1i*2*pi*(1/2*t_tau.^2*B/sampleInterval+fminR*t_tau));
end
% x_direct for now all rows are identical as all direct paths traveled distance * 2 meter 
aoa_direct = 90;


%% STEP 2: Generate multipath signal (x_multipath) at approx 60 deg               %% 
%% x_multipath contains multipath signal at each mic for one chirp duration       %%
%% refer to https://github.com/mikiehan/simulate_fmcw_rx/blob/master/topology.jpg %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%angle of arrival (Same definition as Dina Katabi RF-IDraw)
mpAngles = zeros(1, Nr); % multipath aoa in degree
for i=1:Nr
    mpAngles(i) = atand(distance / (2 * wallDist + spacing * (i-1)));
end
dist3 = sqrt(distance^2 + (2 * wallDist)^2) + [0:Nr-1] * spacing * cosd(mpAngles(1));
dist2 = zeros(1, Nr);
dist1 = zeros(1, Nr);
for i=1:Nr
    dist2(i) = (wallDist + spacing * (i-1))/cosd(mpAngles(i));
    dist1(i) = dist3(i) - dist2(i);
end

att_multipath = zeros(1, Nr); 
for i=1:Nr
     att_multipath(i) = distance^-2 * dist1(i)^-2 * dist2(i)^-2;
end
%att_multipath = ones(1, Nr); %ignore attenuation

distance_multipath = distance + dist3;
tau_multipath = distance_multipath/vs;

t_tau_multipath = zeros(Nr, length(t)); % time stamp for each sample at each microphone for multipath
for i=1:Nr
    t_tau_multipath(i, :) = t - tau_multipath(i);
end

x_multipath = zeros(Nr, length(t));
for i=1:Nr
    %x_multipath(i, :) = att_multipath(i) * sin(2*pi*fc*t_tau_multipath(i, :));
    %x_multipath(i, :) = att_multipath(i) * exp(1i*2*pi*fc*t_tau_multipath(i, :));
    x_multipath(i, :) = att_multipath(i) * exp(1i*2*pi*(1/2*t_tau_multipath(i, :).^2*B/sampleInterval+fminR*t_tau_multipath(i, :)));
end

%% STEP 4-1: Calculate phase difference (mic1 is reference mic)                                       %% 
%% phase_diff contains phase difference of multipath signal among mics                                %%
%% phase_diff2 contains phase difference between direct path and multipath signal at mic1             %%
%% phase_diff + phase_diff2 at row i and sample k is                                                  %%
%% the phase difference at sample k between multipath signal at mic i and direct path signal at mic 1 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


x_multipath_recovered = zeros(Nr, length(t));

phase_diff = zeros(Nr, length(t)); % phase difference among mics for multipath
alpha = B/(2 * sampleInterval);
delay_diff = [0:Nr-1] * spacing * cosd(mpAngles(1)) / vs; % delay difference among mics for multipath (compared to reference mic1)

for i=1:Nr
    phase_diff(i, :) = 2*pi*(-1*fminR * delay_diff(i) + -1 * 2 * alpha * delay_diff(i) * t_tau_multipath(1, :) + alpha * delay_diff(i)^2);
%% debugging to see if multiplying phase_diff indeed can reproduce x_multipath 
%     for k=1:K
%         x_multipath_recovered(i,k) = x_multipath2(i,k) * exp(-1i*phase_diff(i, k)); % should be the same as x_multipath(1,:)
%     end
end

delay_diff2 = tau_multipath(1) - tau; % delay difference between direct path and multipath at mic1 
% phase_diff2 is phase difference between direct path and multipath at reference mic1
phase_diff2 = 2*pi*(-1*fminR * delay_diff2 + -1 * 2 * alpha * delay_diff2 * t_tau(:) + alpha * delay_diff2^2);

x_recovered_mic1 = zeros(1, length(t));
%% debugging to see if multiplying phase_diff2 to direct path signal indeed can reproduce multipath signal received (for mic1)
% for k=1:K
%     x_recovered_mic1(k) = x_multipath_mic1(k) * exp(-1i*phase_diff2(k));%should be the same as x(1,:)
% end

%% STEP 4-2: Based on the phase difference, align all multipath signal  %%
%%  to have the same phase as direct path signal at mic 1               %%
%%  and sum them up                                                     %%
%%  y_MINE_AL contains this sum of aligned signals                      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


w_MINE_MP1 = zeros(Nr, K); % per sample per mic weight
x_recovered = zeros(Nr, K);
for i=1:Nr
    for k=1:K
        w_MINE_MP1(i,k) = exp(-1i*phase_diff(i, k)) * exp(-1i * phase_diff2(k)); % /Nr;
        x_recovered(i, k) = x_multipath(i,k) * w_MINE_MP1(i,k); % undo the phase difference
    end
end
%x_recovered now contains the aligned multipath signal 
% aligned to have the same phase as the direct path signal at mic1

sig_align = [x_direct; x_recovered]; % first Nr rows are direct path sig, the last Nr rows are multipath sig aligned
                              % all Nr * 2 signals have the same phase the direct path signal at mic1
y_MINE_AL = real(sum(sig_align))/(Nr * 2); % add up all aligned signal and normalize


%% STEP 5: Plug in decomposed signal to different beamformin algorithm  %%
%%  such as LCMV, DAS, MVDR and MINE OPT                                %%
%%  where MINE OPT (opt_beam.m) tries to maximize total RSS             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


decomposed_sig = [x_direct; x_multipath];
sig = real(decomposed_sig); % sig contains sum of direct path and multipath
Nr = Nr * 2;

sig = repmat(sig, 1, nChirps);

incidentAz = 90;

y_MINE_AL = repmat(y_MINE_AL, 1, nChirps);

m_xPos = [m_xPos m_xPos];
m_yPos = [m_yPos m_yPos];
m_zPos = [m_zPos m_zPos];
[y_DAS, y_MVDR, y_MVDR2, y_LCMV, y_LP, y_MINE] = beamform(incidentAz, fc, vs, Fs, sig, [], fmaxR, m_xPos, m_yPos, m_zPos, Nr);


figure;
plot(real(y_MINE));

hold on;
plot(real(y_MINE_AL));
plot(real(y_DAS));
plot(real(y_LCMV));
%plot(real(y_MVDR));
%plot(real(y_MVDR2));
plot(real(sig(1, :)));

xlim([0 180]);
legend("MINE no MP", "DAS", "LCMV", "No Beam");


figure;
plot(real(y_MINE));

hold on;

plot(real(y_MINE_AL));
plot(real(sig(1, :)));

xlim([0 180]);
legend("MINE OPT", "MINE ALIGN", "No Beam");


%% STEP 6: De-chirp each beamformed signal and compare the FMCW Profile %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for mic = 1 : Nr
    sig(mic, :) = fftFilter(sig(mic, :),Fs,fminR,fmaxR,50);
end

%figure; spectrogram(Sr_noise(1,:),'yaxis',128,120,128,Fs)

% % de-chirping
[f, profile_mic1] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, sig(1,:), 'No Beam at Mic1');

[f, profile_DAS] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_DAS, 'DAS');
%[f, profile_MVDR2] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_MVDR2, 'MVDR2');
[f, profile_LCMV] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_LCMV, 'LCMV');
%[f, profile_LP] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_LP, 'LP');
[f, profile_MINE] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_MINE, 'MINE OPT');
[f, profile_MINE_AL] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_MINE_AL, 'MINE ALIGN');


figure;
dist = vs*f*sampleInterval*1000/B;
plot(dist, profile_mic1(1,:));
hold on;
plot(dist, profile_DAS(1,:));
%plot(dist, profile_MVDR2(1,:));
plot(dist, profile_LCMV(1,:));
%plot(dist, profile_LP(1,:));
plot(dist, profile_MINE(1,:));
%plot(dist, profile_MINE_MP(1,:));
plot(dist, profile_MINE_AL(1,:));

%plot(dist, profile_FR(1,:));
%plot(dist, profile_PS(1,:));
xlim([0 rtDist*2]);
title ('FMCW Profile (decompose MP)')
xlabel('Distance (m)')
ylabel('Amplitude')
legend("No Beam at Mic1","DAS", "LCMV", "MINE OPT", "MINE ALIGN")

figure;
dist = vs*f*sampleInterval*1000/B;
plot(dist, profile_MINE(1,:));
hold on;
plot(dist, profile_MINE_AL(1,:));
plot(dist, profile_mic1(1,:));
xlim([0 rtDist*2]);
title ('FMCW Profile OPT vs ALIGN (decompose MP)')
xlabel('Distance (m)')
ylabel('Amplitude')
legend("OPT", "ALIGN", "No Beam") %, "PS")

