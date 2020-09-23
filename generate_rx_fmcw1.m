%% ------------------
%% FMCW RX Simulation
%% ------------------
function [Sr_noise, Sr, s_Pos, distance] = generate_rx_fmcw1(fminR, B, Fs, vs, sampleInterval, nChirps, Nr, rxarray, rtDist, xPosWall, mpAngle) 


nChirps = 1;

pos = getElementPosition(rxarray);
m_xPos = pos(1,:); % x-coord
m_yPos = pos(2,:); % y-coord
m_zPos = pos(3,:); % z-coord

s_Pos = [0, rtDist/2, 0]; % source position in meter

xPos = m_xPos(1);
yPos = m_yPos(1);
zPos = m_zPos(1);

distance = sqrt((s_Pos(2) - yPos)^2); % just consider y-axis (one-way direct path distance)


%fminR = 17e3;
%B = 5e3;
%Fs = 48000;
%vs = 340;
%sampleInterval=0.030; % 30 ms
%simLenS =  10;

%Nr = 8; % 8 microphones
%Origin = 4; % 4 meter away approx.
%addNoise = true;
fmaxR = fminR + B;
fc = (fminR + fmaxR)/2;
%lambda = vs/fc; 

Ts=1/Fs;
K=sampleInterval/Ts;

%t=(0:K-1)*Ts;

mpSr = zeros(1, K * (nChirps + 1));

dist2 = sqrt((xPosWall - xPos)^2)/cosd(mpAngle);
delta = tand(mpAngle) * sqrt((xPosWall - xPos)^2);
dist1 = (distance - delta)/sind(mpAngle);
mpAtt = distance^-2 * dist2^-2 * dist1^-2; 
mpTau = (distance + dist1 + dist2)/vs;
    
incidentAngle_mp = [mpAngle; 90];
ebi_mp = steeringVector(m_xPos, m_yPos, m_zPos, fc, vs, incidentAngle_mp(1), incidentAngle_mp(2));
ebi_mp = reshape(ebi_mp, 1, 1, Nr);
ebi_mp_mic1 = ebi_mp(1, 1, 1);

startTic_mp = ceil(mpTau/Ts);
t = (0:K-1)*Ts;
td = t - mpTau;
mpSr0 = mpAtt * cos(2*pi*(1/2*td.^2*B/sampleInterval+fminR*td)) * ebi_mp_mic1;
sig = repmat(mpSr0, 1 , nChirps+1);
mpSr(startTic_mp:end) = sig(1:end-startTic_mp+1);



Sr = zeros(1, K * (nChirps + 1));

att = distance^-4; % attenuation
tau = distance*2/vs; % delay

incidentAngle= [90; 90];
ebi = steeringVector(m_xPos, m_yPos, m_zPos, fc, vs, incidentAngle(1), incidentAngle(2));
ebi = reshape(ebi, 1, 1, Nr);
ebi_mic1 = ebi(1, 1, 1);

startTic = ceil(tau/Ts);
td = t - tau; % new sampling tics with delay
Sr0 = att * cos(2*pi*(1/2*td.^2*B/sampleInterval+fminR*td)) * ebi_mic1; % received Rx signal after tau delay S(t-tau) with attenuation
sig = repmat(Sr0, 1 , nChirps+1);
Sr(startTic:end) = sig(1:end-startTic+1);


Sr = Sr + mpSr; % add multipath signal and direct path signal

Sr_noise = Sr; % For now no noise
% if(addNoise)
%     Sr_noise = Sr + 1e-2 * randn(size(Sr));
% else
%     Sr_noise = Sr;
% end

%cut the initial zeros 
Sr = Sr(startTic:end);
Sr_noise = Sr_noise(startTic:end);

newSize = floor(length(Sr)/K) * K;
Sr = Sr(1:newSize);
Sr_noise = Sr_noise(1:newSize);

end