% https://www.dsprelated.com/freebooks/mdft/Analytic_Signals_Hilbert_Transform.html#fig:sineFD
fminR = 17e3;
B = 5e3;
Fs = 48000;
vs = 340;
sampleInterval=0.030; % 30 ms
%simLenS =  10;

nChirps = 5;

Nr = 4; % 4 microphones
rtDist = 8; % approx. round-trip distance from sound source to microphones
xPosWall = 1; % side wall x position for multipath
mpAngle = 60; % multipath angle in degree
addNoise = true;

fmaxR = fminR + B;
fc = (fminR + fmaxR)/2;

lambda = vs/fc;
spacing = lambda/2;

Ts=1/Fs;
K=sampleInterval/Ts;
simLenS = nChirps * K / Fs;

t=(0:K-1)*Ts;

distance = rtDist/2;
%att = 1; % for now just 1

att = distance^-4; % attenuation

tau =  distance*2/vs;
t_tau = t-tau;

x = zeros(Nr, length(t)); % each row contains direct path signal at mic_row
for i=1:Nr
    %x(i, :) = att * sin(2*pi*fc*t_tau);
    %x(i, :) = att * exp(1i*2*pi*fc*t_tau);
    x(i, :) = att * cos(2*pi*(1/2*t_tau.^2*B/sampleInterval+fminR*t_tau));
end
aoa_direct = 90;

%angle of arrival range 0:180 (Same definition as Dina Katabi RF-IDraw)

aoa_multipath = 60;
x2 = zeros(Nr, length(t)); % each row contains 1 multipath signal at mic_row
d_prime = distance/2/sind(aoa_multipath); % distance at ref mic (last mic)

distance_multipath = spacing*cosd(aoa_multipath)*[0:Nr-1] + d_prime;
att_multipath = ones(1, Nr); % for now just 1
for i=1:Nr
    att_multipath(i) = distance^-2 * distance_multipath(i)^-4;
end
tau_multipath = (distance_multipath * 2 + distance)/vs;

t_tau_multipath = zeros(Nr, length(t));
for i=1:Nr
    t_tau_multipath(i, :) = t - tau_multipath(i);
end

x_multipath = zeros(Nr, length(t));
for i=1:Nr
    %x_multipath(i, :) = att_multipath(i) * sin(2*pi*fc*t_tau_multipath(i, :));
    %x_multipath(i, :) = att_multipath(i) * exp(1i*2*pi*fc*t_tau_multipath(i, :));
    x_multipath(i, :) = att_multipath(i) * cos(2*pi*(1/2*t_tau_multipath(i, :).^2*B/sampleInterval+fminR*t_tau_multipath(i, :)));
end

x_multipath2 = zeros(Nr, length(t));
for i=1:Nr
    %x_multipath2(i, :) = att_multipath(i) * sin(2*pi*fc*t_tau_multipath(1,:)) ...
    %                     * exp(-1i*2*pi*spacing*(i-1)*cosd(aoa_multipath*2)/lambda);
    %x_multipath2(i, :) = att_multipath(i) * exp(1i*2*pi*fc*t_tau_multipath(1,:)) ...
    %                      * exp(-1i*2*pi*spacing*(i-1)*cosd(aoa_multipath)*2/lambda);
    x_multipath2(i, :) = att_multipath(i) * cos(2*pi*(1/2*t_tau_multipath(1, :).^2*B/sampleInterval+fminR*t_tau_multipath(1, :))) ...
                          * exp(-1i*2*pi*spacing*(i-1)*cosd(aoa_multipath)*2/lambda);
end


sig = x + x_multipath; % sig contains sum of direct path and multipath

sig = repmat(sig, 1, nChirps);

incidentAz = 90;
m_xPos = zeros(1, Nr);
m_yPos = spacing * (0:Nr-1);
m_zPos = zeros(1, Nr);

[y_DAS, y_MVDR, y_MVDR2, y_LCMV, y_LP, y_FR] = beamform(incidentAz, fc, vs, Fs, sig, [], fmaxR, m_xPos, m_yPos, m_zPos, Nr);


w = opt_beam(sig, Nr);

y_MINE = w * sig;

figure;
plot(real(y_MINE));

hold on;
plot(real(y_DAS));
plot(real(y_LCMV));
%plot(real(y_MVDR));
plot(real(y_MVDR2));
plot(real(sig(1, :)));

xlim([0 180]);
legend("MINE", "DAS", "LCMV", "MVDR2", "No Beam");
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
[f, profile_MVDR2] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_MVDR2, 'MVDR2');
[f, profile_LCMV] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_LCMV, 'LCMV');
%[f, profile_LP] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_LP, 'LP');
[f, profile_MINE] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, y_MINE, 'MINE');


figure;
dist = vs*f*sampleInterval*1000/B;
plot(dist, profile_mic1(1,:));
hold on;
plot(dist, profile_DAS(1,:));
plot(dist, profile_MVDR2(1,:));
plot(dist, profile_LCMV(1,:));
%plot(dist, profile_LP(1,:));
plot(dist, profile_MINE(1,:));
%plot(dist, profile_FR(1,:));
%plot(dist, profile_PS(1,:));
xlim([0 rtDist*2]);
title ('FMCW Profile')
xlabel('Distance (m)')
ylabel('Amplitude')
legend("No Beam at Mic1","DAS","MVDR2", "LCMV", "MINE") %"FR") %, "PS")
