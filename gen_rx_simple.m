%% ------------------
%% FMCW RX Simulation
%% ------------------
function [Sr_noise, Sr, s_Pos, distance] = gen_rx_simple(fminR, B, Fs, vs, sampleInterval, nChirps, Nr, rxarray, rtDist, xPosWall, mpAngle, addMultipath, addNoise) 


pos = getElementPosition(rxarray);
m_xPos = pos(1,:); % x-coord
m_yPos = pos(2,:); % y-coord
m_zPos = pos(3,:); % z-coord

s_Pos = [0, rtDist/2, 0]; % source position in meter
distance = zeros(1, Nr);

for i=1:Nr
    %X = [s_Pos(1), s_Pos(2); m_xPos(i), m_yPos(i)];
    %distance(i) = pdist(X, 'euclidean');
    distance(i) = sqrt((s_Pos(2) - m_yPos(i))^2); % just consider y-axis (one-way direct path distance)
end

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

mpSr0 = zeros(Nr, K);
mpSr = zeros(Nr, K * (nChirps + 1));


if(addMultipath)

    dist1 = zeros(1, Nr);
    dist2 = zeros(1, Nr);
    mpAtt = zeros(1, Nr);
    mpTau = zeros(1, Nr);
    for i=1:Nr
        dist2(i) = sqrt((xPosWall - m_xPos(i))^2)/cosd(mpAngle);
        delta = tand(mpAngle) * sqrt((xPosWall - m_xPos(i))^2);
        dist1(i) = (distance(i) - delta)/sind(mpAngle);
        mpAtt(i) = distance(i)^-2 * dist2(i)^-2 * dist1(i)^-2; 
        mpTau(i) = (distance(i) + dist1(i) + dist2(i))/vs;
    end
    
    
    incidentAngle = [mpAngle*-1; 90];
    ebi = steeringVector(m_xPos, m_yPos, m_zPos, fc, vs, incidentAngle(1), incidentAngle(2));
    ebi = reshape(ebi, 1, 1, Nr);

    for i=1:Nr
        startTic = ceil(mpTau(i)/Ts);
        t = (0:K-1)*Ts;
        td = t - mpTau(i);
        %td = t;
        mpSr0(i, :) = mpAtt(i) * cos(2*pi*(1/2*td.^2*B/sampleInterval+fminR*td)) * ebi(1,1, i);
        sig = repmat(mpSr0(i, :), 1 , nChirps+1);
        mpSr(i, startTic:end) = sig(1:end-startTic+1);
    end
    
end


att = zeros(1, Nr);
tau = zeros(1, Nr);
Sr0 = zeros(Nr, K);
Sr = zeros(Nr, K * (nChirps + 1));
%Sr_noise = zeros(Nr, K * nChirps);

for i=1:Nr
    att(i) = distance(i)^-4; % attenuation
    tau(i) = distance(i)*2/vs; % delay
end


for i=1:Nr
    startTic = ceil(tau(i)/Ts);
    t=(0:K-1)*Ts;
    td = t - tau(i); % new sampling tics with delay
    %td = t;
    Sr0(i, :) = att(i) * cos(2*pi*(1/2*td.^2*B/sampleInterval+fminR*td)); % received Rx signal after tau delay S(t-tau) with attenuation
    sig = repmat(Sr0(i, :), 1 , nChirps+1);
    Sr(i, startTic:end) = sig(1:end-startTic+1);
end

if(addMultipath)
    Sr = Sr + mpSr;
end

if(addNoise)
    Sr_noise = Sr + 1e-2 * randn(size(Sr));
else
    Sr_noise = Sr;
end

% cutoff the first chirp 
Sr = Sr(:, K+1:end);
Sr_noise = Sr_noise(:, K+1:end);

end