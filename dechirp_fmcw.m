function [f, profile] = dechirp_fmcw(rtDist, Fs, fminR, B, vs, sampleInterval, nChirps, Sr, algo_name)
%fminR = 17e3;
%B = 5e3;
%Fs = 48000;
%vs = 340;
%sampleInterval=0.030; % 30 ms

Ts=1/Fs;
K=sampleInterval/Ts;

t=(0:K-1)*Ts;

St0 = cos(2*pi*(1/2*t.^2*B/sampleInterval+fminR*t)); % original Tx signal
St = repmat(St0, 1 , nChirps);

%synchronizing to remove the intial noise
% maxIndex=syncFMCWSymbol(Sr,0,St0,length(St0));
% 
% Sr = Sr(maxIndex:length(Sr));
% St = St(1:length(Sr));

Sm = St .* Sr; % mixing St and Sr
% de-chirp
figure;
profile = [];
for c=1:nChirps-1
    [f,fft1] = plotFFT(Sm(1 + (c-1) * K : c*K),Fs);
    dist = vs*f*sampleInterval*1000/B;
    plot(dist, fft1); % should peak at D
    hold on;
    xlim([0 rtDist*2]);
    %ylim([0 0.01]);
    profile = [profile; fft1];
end
title(algo_name)
xlabel('Distance (m)')
ylabel('Amplitude')
end