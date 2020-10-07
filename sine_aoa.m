fminR = 17e3;
B = 5e3;
Fs = 5000;
vs = 340;
sampleInterval=0.050; % 50 ms
simLenS =  10;

fmaxR = fminR + B;
fc = (fminR + fmaxR)/2;

fc = 50;
lambda = vs/fc; 
spacing = lambda/2;

Ts=1/Fs;
K=sampleInterval/Ts;

t=(0:K-1)*Ts;

wave = sin(2*pi*fc*t);

wave2 = exp(1i*2*pi*fc*t);

figure;
plot(t, real(wave2));
hold on;
plot(t, imag(wave2));
legend("real", "imag");


wave3 = wave2 * exp(-1i*2*pi*0.125); % negative
wave4 = exp(1i*2*pi*fc*(t-0.0025)); % longer 

figure;
plot(t, real(wave2));
hold on;
%plot(t, imag(wave2));
plot(t, real(wave3));
plot(t, real(wave4));
%plot(t, imag(wave3));
legend("real", "real3-", "real4ph");
%legend("real", "imag", "real3", "imag3");


figure;
plot(t, imag(wave2));
hold on;
%plot(t, imag(wave2));
plot(t, imag(wave3));
plot(t, imag(wave4));
%plot(t, imag(wave3));
legend("imag", "imag3-", "imag4ph");
%legend("real", "imag", "real3", "imag3");
