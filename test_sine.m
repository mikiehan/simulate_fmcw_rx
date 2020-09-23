% https://www.dsprelated.com/freebooks/mdft/Analytic_Signals_Hilbert_Transform.html#fig:sineFD
fminR = 17e3;
B = 5e3;
Fs = 48000;
vs = 340;
sampleInterval=0.030; % 30 ms
simLenS =  10;

Nr = 2; % 4 microphones
rtDist = 8; % approx. round-trip distance from sound source to microphones
xPosWall = 1; % side wall x position for multipath
mpAngle = 60; % multipath angle in degree
addNoise = true;

fmaxR = fminR + B;
fc = (fminR + fmaxR)/2;
fc = 500;
lambda = vs/fc; 
spacing = lambda/2;

Ts=1/Fs;
K=sampleInterval/Ts;

t=(0:K-1)*Ts;

distance = rtDist/2;
att = 1; % for now just 1

att = distance^-4; % attenuation

tau =  distance*2/vs;
t_tau = t-tau;

x = zeros(Nr, length(t)); % each row contains direct path signal at mic_row
for i=1:Nr
    x(i, :) = att * sin(2*pi*fc*t_tau);
    %x(i, :) = att * exp(1i*2*pi*fc*t_tau);
end

%angle of arrival range 0:180 (Same definition as Dina Katabi RF-IDraw)
aoa_direct = 90; 
aoa_multipath = 60;
x2 = zeros(Nr, length(t)); % each row contains 1 multipath signal at mic_row
d_prime = distance/2/sind(aoa_multipath); % distance at ref mic (last mic)

distance_multipath = spacing*cosd(aoa_multipath)*[0:Nr-1] + d_prime;
att_multipath = ones(1, Nr); % for now just 1
for i=1:Nr
    att_multipath(i) = distance^-2 * distance_multipath(i)^-4; 
end
tau_multipath = (distance_multipath * 2 + distance)/vs;

t_tau_multipath = zeros(Nr, length(t));
for i=1:Nr
    t_tau_multipath(i, :) = t - tau_multipath(i);
end

%x_multipath = zeros(Nr, length(t));
%for i=1:Nr
%    x_multipath(i, :) = att_multipath(i) * sin(2*pi*fc*t_tau_multipath(i, :));
%    %x_multipath(i, :) = att_multipath(i) * exp(1i*2*pi*fc*t_tau_multipath(i, :));
%end

x_multipath2 = zeros(Nr, length(t));
for i=1:Nr
    x_multipath2(i, :) = att_multipath(i) * sin(2*pi*fc*t_tau_multipath(1,:)) ...
                         * exp(1i*2*pi*spacing*cosd(aoa_multipath)/lambda*(i-1));
    %x_multipath2(i, :) = att_multipath(i) * exp(1i*2*pi*fc*t_tau_multipath(1,:)) ...
    %                     * exp(1i*2*pi*spacing*cosd(aoa_multipath)/lambda*(i-1));
end 


sig = x + x_multipath2; % sig contains sum of direct path and multipath 

opt_beam(sig, Nr)
