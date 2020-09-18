%% ------------------
%% FMCW RX Simulation
%% ------------------
function [m_xPos, m_yPos, m_zPos, rxarray]  = generate_rx_tx_positions(Nr, radius, fmaxR)

%Nr = 8; % 8 microphones
%Origin = 4; % 4 meter away approx.
%fmaxR = fminR + B;
%radius = 0.05;  % array radius (m)

microphone = phased.OmnidirectionalMicrophoneElement('FrequencyRange',[20 fmaxR]);
rxarray = phased.UCA('NumElements', Nr ,'Radius', radius, "Element", microphone);
%rxarray = phased.ULA('NumElements',Nr,'ElementSpacing',radius, "Element", microphone);
% pos0 = getElementPosition(rxarray0); 
% xPos = pos0(2,:);
% yPos = pos0(1,:);
% zPos = pos0(3,:);
% 
% rxarray = phased.ConformalArray(...
%     'ElementPosition',[xPos;yPos;zPos],...
%     'Element', microphone);
%     %'ElementNormal',[ones(1,Nr)*90;zeros(1,Nr)],...

figure;
viewArray(rxarray, 'ShowNormals', true);

pos = getElementPosition(rxarray);
%normvec = getElementNormal(rxarray);
m_xPos = pos(1,:); % x-coord
m_yPos = pos(2,:); % y-coord
m_zPos = pos(3,:); % z-coord


end
