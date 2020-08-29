%% -------------------
%% Radar Signal Simulation
%% -------------------
fminR = 18e3;
B = 2e3;
Fs = 48000;
vs = 340;
sampleInterval=0.030; % 30 ms
simLenS =  10;

Ts=1/Fs;
K=sampleInterval/Ts;


t=(0:K-1)*Ts;
Fi=(0:K-1)*(B/K)+fminR;

yR0=cos(2*pi*(1/2*t.^2*B/sampleInterval+fminR*t));


nChirps = floor(simLenS * Fs / K);
simLenS = nChirps * K / Fs;

tx = repmat(yR0, 1, nChirps);

tx = fftFilter(tx,Fs,fminR,fminR+B,50);
figure; spectrogram(tx,'yaxis',128,120,128,Fs)%plotSignal(TX,Fs);

fi = repmat(Fi, 1, nChirps); % instantaneous FMCW freq
rx = zeros(size(tx)); 
dist = 2; % round trip distance
lag = floor(dist / vs * Fs);

for ti = 1:length(rx)
    if ti > lag
        %mag = (1/(dist/2))^4;
        mag = 1; % for now no attenuation due to distance
         % simulated rx (adding attenuation, phase shift and delay) 
        rx(ti) = mag*exp(-1i*2*pi*fi(ti-lag)*dist/vs)*tx(ti - lag);
    end
end

% FMCW de-chirping code start 
y = fftFilter(rx,Fs,fminR,fminR+B,500);

% TX signal
TX = genChirp(Fs,sampleInterval,0,fminR,B,length(y)).';

%synchronizing
yR0 = genChirp(Fs,sampleInterval,0,fminR,B,1);

maxIndex=syncFMCWSymbol(y,0,yR0,length(yR0));
disp(maxIndex);
y = y(maxIndex:length(y));

figure; spectrogram(y,'yaxis',128,120,128,Fs)

% windowing RX signal
w=0.5-0.5*cos(2*pi/K*(0:length(y)-1));
y = y.*w;
figure; spectrogram(y(1:Fs*1),'yaxis',128,120,128,Fs)


% de-chirping
results = [];
Ne = floor(length(y)/K)-1;
figure; % figure should show clear peak at distance 
for i = 1:Ne % for each chirp
    y0 = y(1+(i-1)*K:i*K+1);
    x0 = TX(1+(i-1)*K:i*K+1).';
    prod = y0.*x0;
    [f,fft1] = plotFFT(prod,Fs); 
    dist = vs*f*sampleInterval*1000/B;

    plot(dist, fft1);
    hold on;
    xlim([0 3]); 
    xlabel('distance (m)')
    results = [results; fft1];
end
