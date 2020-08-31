%% ------------------
%% FMCW RX Simulation
%% ------------------
function [Sr_noise, Sr] = generate_rx_fmcw(fminR, B, Fs, vs, sampleInterval, nChirps, Nr, distance, addNoise) 

%fminR = 17e3;
%B = 5e3;
%Fs = 48000;
%vs = 340;
%sampleInterval=0.030; % 30 ms
%simLenS =  10;

%Nr = 8; % 8 microphones
%Origin = 4; % 4 meter away approx.
%addNoise = true;

Ts=1/Fs;
K=sampleInterval/Ts;

t=(0:K-1)*Ts;

att = zeros(1, Nr);
tau = zeros(1, Nr);
Sr0 = zeros(Nr, K);
Sr = zeros(Nr, K * nChirps);
Sr_noise = zeros(Nr, K * nChirps);

for i=1:Nr
    att(i) = (distance(i)/2)^-4; % attenuation
    tau(i) = distance(i)/vs; % delay
    td = t - tau(i); % new time tics with delay
    Sr0(i, :) = att(i) * cos(2*pi*(1/2*td.^2*B/sampleInterval+fminR*td)); % received Rx signal after tau delay S(t-tau) with attenuation
    Sr(i, :) = repmat(Sr0(i, :), 1 , nChirps);
end

if(addNoise)
    Sr_noise = Sr + 1e-2 * randn(size(Sr));
else
    Sr_noise = Sr;
end

end