% https://www.dsprelated.com/freebooks/mdft/Analytic_Signals_Hilbert_Transform.html#fig:sineFD
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
addNoise = false;
SNR_db = 5;

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


att = distance^-4; % attenuation
%att = 1; %ignore attenuation
tau =  distance*2/vs;
t_tau = t-tau;

x = zeros(Nr, length(t)); % each row contains direct path signal at mic_row
for i=1:Nr
    %x(i, :) = att * sin(2*pi*fc*t_tau);
    %x(i, :) = att * exp(1i*2*pi*fc*t_tau);
    x(i, :) = att * exp(1i*2*pi*(1/2*t_tau.^2*B/sampleInterval+fminR*t_tau));
end
aoa_direct = 90;

%angle of arrival range 0:180 (Same definition as Dina Katabi RF-IDraw)
mpAngles = zeros(1, Nr); % multipath aoa in degree
for i=1:Nr
    mpAngles(i) = atand(distance / (2 * wallDist + spacing * (i-1)));
end
%dist3 = sqrt(distance^2 + (2 * wallDist + spacing * [0:Nr-1]).^2); % distance of multipath coming back from object to mic
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

t_tau_multipath = zeros(Nr, length(t));
for i=1:Nr
    t_tau_multipath(i, :) = t - tau_multipath(i);
end

x_multipath = zeros(Nr, length(t));
for i=1:Nr
    %x_multipath(i, :) = att_multipath(i) * sin(2*pi*fc*t_tau_multipath(i, :));
    %x_multipath(i, :) = att_multipath(i) * exp(1i*2*pi*fc*t_tau_multipath(i, :));
    x_multipath(i, :) = att_multipath(i) * exp(1i*2*pi*(1/2*t_tau_multipath(i, :).^2*B/sampleInterval+fminR*t_tau_multipath(i, :)));
end

x_multipath2 = zeros(Nr, length(t));
x_multipath_recovered = zeros(Nr, length(t));

phase_diff = zeros(Nr, length(t)); % phase difference among mics for multipath
alpha = B/(2 * sampleInterval);
delay_diff = [0:Nr-1] * spacing * cosd(mpAngles(1)) / vs; % delay difference among mics (compared to reference mic1)

for i=1:Nr
    phase_diff(i, :) = 2*pi*(-1*fminR * delay_diff(i) + -1 * 2 * alpha * delay_diff(i) * t_tau_multipath(1, :) + alpha * delay_diff(i)^2);
%% debugging to see if multiplying phase_diff indeed can reproduce x_multipath 
%     for k=1:K
%         x_multipath2(i, k) = x_multipath(1,k) * exp(1i*phase_diff(i, k));
%         x_multipath_recovered(i,k) = x_multipath2(i,k) * exp(-1i*phase_diff(i, k)); % should be the same as x_multipath(1,:)
%     end
end

delay_diff2 = tau_multipath(1) - tau; % delay difference between direct path and multipath at mic1 
% phase_diff2 is phase difference between direct path and multipath at reference mic1
phase_diff2 = 2*pi*(-1*fminR * delay_diff2 + -1 * 2 * alpha * delay_diff2 * t_tau(:) + alpha * delay_diff2^2);

x_multipath_mic1 = zeros(1, length(t));
x_recovered_mic1 = zeros(1, length(t));

%% debugging to see if multiplying phase_diff2 to direct path signal indeed can reproduce multipath signal received (for mic1)
% for k=1:K
%     x_multipath_mic1(k) = x(1,k) * exp(1i * phase_diff2(k)); 
%     x_recovered_mic1(k) = x_multipath_mic1(k) * exp(-1i*phase_diff2(k));
% end

w_MINE_MP1 = zeros(Nr, K); % per sample per mic weight
x_recovered_all = zeros(Nr, K);
x_recovered_all_approx = zeros(Nr, K);
for i=1:Nr
    for k=1:K
        w_MINE_MP1(i,k) = exp(-1i*phase_diff(i, k)) * exp(-1i * phase_diff2(k)); % /Nr;
        x_recovered_all(i, k) = x_multipath(i,k) * w_MINE_MP1(i,k);
    end
end

figure;
plot(real(w_MINE_MP1(1,:)));

figure;
plot(real(w_MINE_MP1(2,:)));

% figure;
% plot(imag(w_MINE_MP(1,:)));
% 
% figure;
% plot(imag(w_MINE_MP(2,:)));

% decomposed_sig = [x;x_multipath]; % direct path signal x and multipath signal x_recovered
% w_MINE_ALL = opt_beam(decomposed_sig, Nr*2);
% w_MINE_ALL = w_MINE_ALL/sum(abs(w_MINE_ALL)); %normalize weight
% y_MINE_ALL = w_MINE_ALL * decomposed_sig;
% 
% w_MINE_MP2 = opt_beam(x_multipath, Nr);
% w_MINE_MP2 = w_MINE_MP2/sum(abs(w_MINE_MP2)); %normalize weight
% y_MINE_MP = w_MINE_MP2 * x_recovered_all; 

y_MINE_MP = real(sum(x_recovered_all))/Nr; 
y_MINE_ALL = real(sum(x + x_recovered_all))/(2*Nr);

sig = real(x + x_multipath); % sig contains sum of direct path and multipath


sig = repmat(sig, 1, nChirps);
noise = zeros(size(sig));
if(addNoise) 
    P0 = mean(sig.^2);
    SNR = 10^(SNR_db/10);
    noise = sqrt(P0/SNR).*randn(size(sig));
    sig = sig + noise;
end

incidentAz = 90;

%sig2 = [x ; x_multipath];
%sig2 = repmat(sig2, 1, nChirps);

y_MINE_MP = repmat(y_MINE_MP, 1, nChirps);
y_MINE_ALL = repmat(y_MINE_ALL, 1, nChirps);

% w_MINE2 = opt_beam(sig2, Nr * 2);
% w_MINE2 = w_MINE2/sum(abs(w_MINE2)); % normalize weight
% y_MINE2 = w_MINE2 * sig2;
% 
% sig3 = repmat(x_multipath, 1, nChirps);
% w_MINE3 = opt_beam(sig3, Nr);
% w_MINE3 = w_MINE3/sum(abs(w_MINE3)); % normalize weight
% y_MINE3 = w_MINE3 * sig3;
% 
% sig4 = repmat(x, 1, nChirps);
% w_MINE4 = opt_beam(sig4, Nr);
% w_MINE4 = w_MINE4/sum(abs(w_MINE4)); % normalize weight
% y_MINE4 = w_MINE4 * sig4;
% 
% sig5 = [y_MINE4; y_MINE3];
% w_MINE5 = opt_beam(sig5, 2);
% w_MINE5 = w_MINE5/sum(abs(w_MINE5)); % normalize weight
% y_MINE5 = w_MINE5 * sig5;
% 

[y_DAS, y_MVDR, y_MVDR2, y_LCMV, y_LP, y_MINE] = beamform(incidentAz, fc, vs, Fs, sig, [], fmaxR, m_xPos, m_yPos, m_zPos, Nr);


figure;
plot(real(y_MINE));

hold on;
plot(real(y_MINE_MP));
plot(real(y_MINE_ALL));
plot(real(y_DAS));
plot(real(y_LCMV));
%plot(real(y_MVDR));
%plot(real(y_MVDR2));
plot(real(sig(1, :)));

xlim([0 180]);
legend("MINE no MP", "MINE MP", "DAS", "LCMV", "No Beam");


figure;
plot(real(y_MINE));

hold on;
plot(real(y_MINE_MP));

plot(real(y_MINE_ALL));
%plot(real(y_MVDR));
%plot(real(y_MVDR2));
plot(real(sig(1, :)));

xlim([0 180]);
legend("MINE", "MINE MP", "MINE ALL", "No Beam");


%legend("MINE", "DAS", "MVDR2", "No Beam");

% figure;
% plot(real(y_MVDR2));
% hold on;
% plot(real(sig(1, :)));
% xlim([0 180]);
% legend("MVDR2", "No Beam");

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
[f, profile_MINE] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_MINE, 'MINE');
[f, profile_MINE_MP] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_MINE_MP, 'MINE MP');
[f, profile_MINE_ALL] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_MINE_ALL, 'MINE ALL');


figure;
dist = vs*f*sampleInterval*1000/B;
plot(dist, profile_mic1(1,:));
hold on;
plot(dist, profile_DAS(1,:));
%plot(dist, profile_MVDR2(1,:));
plot(dist, profile_LCMV(1,:));
%plot(dist, profile_LP(1,:));
plot(dist, profile_MINE(1,:));
plot(dist, profile_MINE_MP(1,:));
plot(dist, profile_MINE_ALL(1,:));

%plot(dist, profile_FR(1,:));
%plot(dist, profile_PS(1,:));
xlim([0 rtDist*2]);
title ('FMCW Profile')
xlabel('Distance (m)')
ylabel('Amplitude')
legend("No Beam at Mic1","DAS", "LCMV", "MINE", "MINE MP", "MINE ALL")

figure;
dist = vs*f*sampleInterval*1000/B;
plot(dist, profile_MINE(1,:));
hold on;
plot(dist, profile_MINE_MP(1,:));
plot(dist, profile_MINE_ALL(1,:));
plot(dist, profile_mic1(1,:));
xlim([0 rtDist*2]);
title ('FMCW Profile Considering No MP, only MP, both')
xlabel('Distance (m)')
ylabel('Amplitude')
legend("No MP", "Only MP", "BOTH", "No Beam") %, "PS")

