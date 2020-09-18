fs = 1000; % Sampling frequency (samples per second)
dt = 1/fs; % seconds per sample
StopTime = 0.25; % seconds
t = (0:dt:StopTime)'; % seconds
F = 60; % Sine wave frequency (hertz)
data = sin(2*pi*F*t);
t_tau = t - 0.025;
data_tau = sin(2*pi*F*t_tau);

figure;
plot(t,data)
hold on;
plot(t_tau, data_tau);

%%For one cycle get time period
T = 1/F ;
% time step for one time period
tt = 0:dt:T+dt ;
d = sin(2*pi*F*tt) ;
tt_tau = tt - 0.025; % 0.0025 delay
d_tau = sin(2*pi*F*tt_tau);
figure;
%plot(tt,d);
hold on;
plot(tt_tau, d_tau);

x = 0.125;
%d_tau_translate = d_tau * exp(1i*2*pi*x);

d_tau_translate = zeros(1,length(d_tau));
for i=1:length(d_tau)
    d_tau_translate(i) = d_tau(i) * (cosd(60) - cos(2*pi*F*tt_tau(i)) * sind(60) / d_tau(i));
end

plot(tt_tau, d_tau_translate);

d_tau_translate2 = zeros(1,length(d_tau));
for i=1:length(d_tau)
    d_tau_translate2(i) = d_tau(i) * exp(1i*pi/3);
end

plot(tt_tau, d_tau_translate2);
legend("org", "t60", "t60-exp");

%d_tau_translate2 = d_tau * (cosd(90) - cos(2*pi*F*tt_tau) * sind(90) / d_tau);
%plot(tt_tau, d_tau_translate2);
%legend("org","trans 60" , "trans 90");
