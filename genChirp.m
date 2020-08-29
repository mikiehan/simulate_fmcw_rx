function [ signal ] = genChirp( Fs,sampleInterval,offInterval,fmin,B,len )

Ts=1/Fs;
K=sampleInterval/Ts;
K=round(K);

K2 = offInterval/Ts;
K2=round(K2);
z = zeros(1,K2);

t=(0:K-1)*Ts;

% generate a single FMCW signals
%yR0=exp(-1i*2*pi*(1/2*t.^2*B/sampleInterval+fmin*t));
yR0=cos(2*pi*(1/2*t.^2*B/sampleInterval+fmin*t));


%FMCWSymbolRef=exp(-1i*2*pi*(1/2*te.^2*B/FMCWLength+fmin*te)).';


% repeat the single FMCW to generate 
sig = [yR0,z];

numReps = ceil(len/length(sig));
signal = repmat(sig,1,numReps);
%figure; spectrogram(signal,'yaxis',128,120,128,Fs)%plotSignal(TX,Fs);
end

