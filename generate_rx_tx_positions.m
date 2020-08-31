%% ------------------
%% FMCW RX Simulation
%% ------------------
function [s_Pos, m_xPos, m_yPos, m_zPos,rxarray,distance]  = generate_rx_tx_positions(Nr, radius, sOrigin, fmaxR)

%Nr = 8; % 8 microphones
%Origin = 4; % 4 meter away approx.
%fmaxR = fminR + B;
%radius = 0.05;  % array radius (m)

microphone = phased.OmnidirectionalMicrophoneElement('FrequencyRange',[20 fmaxR]);
rxarray = phased.UCA('NumElements', Nr ,'Radius', radius, "Element", microphone);
pos = getElementPosition(rxarray);
m_xPos = pos(1,:); % x-coord
m_yPos = pos(2,:); % y-coord
m_zPos = pos(3,:); % z-coord

s_Pos = [0, sOrigin, 0]; % source position in meter
distance = zeros(1, Nr);

for i=1:Nr
    X = [s_Pos(1), s_Pos(2); m_xPos(i), m_yPos(i)];
    distance(i) = pdist(X, 'euclidean');
end

end