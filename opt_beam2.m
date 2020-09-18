function [x, w] = opt_beam2(Sr) % for now num mics are hard-coded to 8

fminR = 17e3;
B = 5e3;
Fs = 48000;
%vs = 340;
sampleInterval=0.030; % 30 ms

Ts=1/Fs;
K=sampleInterval/Ts;
t = (0:K-1)*Ts;

lb = [];
ub = [];

sig = cos(2*pi*(1/2*t.^2*B/sampleInterval+fminR*t));
objective = @(x) func_obj(x, sig);
nonlincon = @(x) constraint(x, sig, Sr);

x0 = [1 1 1];

x = fmincon(objective,x0,[],[],[],[],lb,ub,nonlincon);

disp(['x1 = ' num2str(x(1))])
disp(['x2 = ' num2str(x(2))])
disp(['x3 = ' num2str(x(3))])

%w = exp(1i*x); %weights

% x1 is coefficient for direct path signal
% x2 is coefficient for multi path signal
% x3 is weight for x2 

w = exp(1i*x(3)); % weight

% Maximize RSS (sum of each sample's absolute value)
% function obj = func_obj(x,r)
%    obj = sum(abs(exp(1i*x(1))*r(1, :)+ exp(1i*x(2))*r(2, :) + exp(1i*x(3))*r(3, :) ...
%          + exp(1i*x(4))*r(4, :)+ exp(1i*x(5))*r(5, :) + exp(1i*x(6))*r(6, :) ...
%          + exp(1i*x(7))*r(7, :)+ exp(1i*x(8))*r(8, :)))*-1;
% end

    function obj = func_obj(x, sig)
        obj = sum(abs(x(1) * sig + x(2) * sig * exp(1i*x(3)))); 
    end


end

