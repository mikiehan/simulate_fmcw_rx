function [y_DAS, y_MVDR, y_LCMV, y_LP, y_FR] = beamform(incidentAz, fc, c, fs, Sr, rxarray, fmaxR, m_xPos, m_yPos, m_zPos, num_mics)

x = Sr'; 
x = x(:,1:num_mics);

% parameters
%c = 340;        % sound speed
%fs = 48000;     % sampling rate
%duration = 5;   % duration time
%Nr = 8;         % # of antenna element of Rx
%Ns = 1;         % # source signal
%incidentAz = 90; % incident angle (from where sound source is coming)

m_xPos_used = m_xPos(1:num_mics);
m_yPos_used = m_yPos(1:num_mics);
m_zPos_used = m_zPos(1:num_mics);


thetaScanningAngles = 0:1:180; %0:1:360;
phiScanningAngles = 90; %0:1:180; % 90
escan = steeringVector(m_xPos_used, m_yPos_used, m_zPos_used, fc, c, thetaScanningAngles, phiScanningAngles);

%%
incidentAngle = [incidentAz; 90];
% plot beamform patterns
ebi = steeringVector(m_xPos_used, m_yPos_used, m_zPos_used, fc, c, incidentAngle(1), incidentAngle(2));
ebi = reshape(ebi, 1, 1, num_mics);


% beamformer = phased.PhaseShiftBeamformer('SensorArray',rxarray,...
%     'OperatingFrequency',fc,'PropagationSpeed',c,...
%     'Direction',[incidentAz; 0],'WeightsOutputPort',true);
% [y_PS,w_PS] = beamformer(x);

normvec = getElementNormal(rxarray);

microphone = phased.OmnidirectionalMicrophoneElement('FrequencyRange',[20 fmaxR]);
rxarray_trimmed = phased.ConformalArray(...
    'ElementPosition',[m_xPos_used;m_yPos_used;m_zPos_used],...
    'Element', microphone, ...
    'ElementNormal',normvec(:,num_mics));

beamformer = phased.FrostBeamformer('SensorArray',rxarray_trimmed,...
     'PropagationSpeed',c,'SampleRate', fs, ...
     'Direction',[incidentAz; 0],'WeightsOutputPort',true);
 [y_FR,w_FR] = beamformer(x);

w_DAS = weightingVectorDAS(ebi);
w_MVDR = weightingVectorMVDR(x.', ebi);
w_LCMV = weightingVectorLCMV(x.');
w_LP = weightingVectorLP(x.');
AF_DAS = arrayFactor(w_DAS, escan);
AF_MVDR = arrayFactor(w_MVDR, escan);
AF_LCMV = arrayFactor(w_LCMV, escan);
AF_LP = arrayFactor(w_LP, escan);
AF_FR = arrayFactor(w_FR, escan);
% AF_PS = arrayFactor(w_PS, escan);

figure;
polarplot(deg2rad(thetaScanningAngles),AF_DAS')
hold on
polarplot(deg2rad(thetaScanningAngles),AF_MVDR')
polarplot(deg2rad(thetaScanningAngles),AF_LCMV')
polarplot(deg2rad(thetaScanningAngles),AF_LP')
polarplot(deg2rad(thetaScanningAngles),AF_FR')
%polarplot(deg2rad(thetaScanningAngles),AF_PS')
hold off

% polarplot(deg2rad(phiScanningAngles),AF_DAS')
% hold on
% polarplot(deg2rad(phiScanningAngles),AF_MVDR')
% polarplot(deg2rad(phiScanningAngles),AF_LCMV')
% polarplot(deg2rad(phiScanningAngles),AF_LP')
% % polarplot(deg2rad(phiScanningAngles),AF_PS')
% hold off

%legend(["DAS", "MVDR"]);
legend(["DAS", "MVDR", "LCMV", "LP", "FR"]); %, "PS"])
%%
% Beamforming
y_DAS = (x*squeeze(w_DAS))';
y_MVDR = (x*squeeze(w_MVDR))';
y_LCMV = (x*squeeze(w_LCMV))';
y_LP = (x*squeeze(w_LP))';
y_FR = y_FR';
%y_PS = y_PS';
% figure(4);
% plot(t,real(x(:,1)), t,real(x_raw(:,1)));
% xlim([0 0.001]);
% legend('Raw','Raw no noise')

% t = 0:1/Fs:duration;
% t = t(1:end-1)';
% 
% figure(4);
% plot(t,real(x(:,2)),'b:',t,real(y_DAS),'g', t, real(y_MVDR),'c', t, real(y_LCMV), 'y', t, real(y_LP), 'r');
% xlim([0 0.005]);
% xlabel('Time (s)')
% ylabel('Amplitude')
% legend('Raw','DAS','MVDR', "LCMV", "LP")

%%
% calculate array gain
% snr_raw = snr(real(x(:,2)), noise(:,2))
% snr_DAS = snr(real(y_DAS), noise(:,2))
% snr_MVDR = snr(real(y_MVDR), noise(:,2))
% snr_LCMV = snr(real(y_LCMV), noise(:,2))
% snr_LP = snr(real(y_LP), noise(:,2))
% snr_PS = snr(real(y_PS), noise(:,2))
% 
% disp(snr_DAS - snr_raw);
% disp(snr_LCMV - snr_raw);
% disp(snr_LP - snr_raw);

% figure(5);
% plot(t,real(x(:,2)),'b:',t,real(y_DAS), 'g', t, real(y_LCMV), 'c', t, real(y_LP), 'r');
% xlim([0 0.001]);
% xlabel('Time (s)')
% ylabel('Amplitude')
% legend('Raw','DAS', "LCMV", "LP")
% 
% normalized_yDAS = real(y_DAS)./max(abs(real(y_DAS)));
% normalized_yLCMV = real(y_LCMV)./max(abs(real(y_LCMV)));
% normalized_yLP = real(y_LP)./max(abs(real(y_LP)));

% figure; spectrogram(y_DAS,'yaxis',128,120,128,fs);
% %audiowrite(strcat(dir , prefix , '_channel_8_das.wav'), real(y_DAS), fs);
% audiowrite(strcat(dir , prefix , '_channel_8_das.wav'), normalized_yDAS, fs);
% 
% %audiowrite(strcat(dir , prefix , '_channel_8_lcmv.wav'), real(y_LCMV), fs);
% audiowrite(strcat(dir , prefix , '_channel_8_lcmv.wav'), normalized_yLCMV, fs);
% %audiowrite(strcat(dir , prefix , '_channel_8_lp.wav'), real(y_LP), fs);
% audiowrite(strcat(dir , prefix , '_channel_8_lp.wav'), normalized_yLP, fs);
% 
% audiowrite(strcat(dir , prefix , '_channel_0.wav'), rx_raw(:,1), fs);
