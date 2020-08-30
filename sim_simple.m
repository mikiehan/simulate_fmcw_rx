%% ------------------
%% FMCW RX Simulation
%% ------------------
fminR = 20e3;
B = 5e3;
Fs = 48000;
vs = 340;
sampleInterval=0.030; % 30 ms
simLenS =  10;

Ts=1/Fs;
K=sampleInterval/Ts;

t=(0:K-1)*Ts;
nChirps = floor(simLenS * Fs / K);
simLenS = nChirps * K / Fs;


St = cos(2*pi*(1/2*t.^2*B/sampleInterval+fminR*t)); % original Tx signal

D = 4;
tau = D/vs; % delay
td = t - tau;
att = (D/2)^-4; % half of distance 
att = 1;
Sr = att * cos(2*pi*(1/2*td.^2*B/sampleInterval+fminR*td)); % received Rx signal after tau delay S(t-tau) with attenuation

Sm = St .* Sr; % mixing St and Sr 
% Should be approx cos(2*pi*fminR*(D/vs) + 2*pi*B*t*(D/(vs*sampleInterval)) - pi*B*(D*D/(vs*vs*sampleInterval)));

% de-chirp
[f,fft1] = plotFFT(Sm,Fs); 
figure; 
dist = vs*f*sampleInterval*1000/B;
plot(dist, fft1); % should peak at D 
max(fft1)

