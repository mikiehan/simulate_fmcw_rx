function w = weightingVectorMVDR(inputSignal, e)
%weightingVectorMinimumVariance - minimum variance optimal weights
%
%Calculates the optimal minimum variance weights for a specific array with
%a specific input signal and specific scanning angle
%
%w = weightingVectorMinimumVariance(inputSignal, e)
%
%IN
%inputSignal         - PxL vector of inputsignals consisting of L samples
%e                   - MxNxP steering vector of P antennas

%
%OUT
%w                   - MxNxP vector of complex optimal weights
%

%M # of theta steering angles, N # of phi steering angless, P number of mics
[M, N, P] = size(e);

%Calculate correlation matrix with diagonal loading
R = inputSignal*inputSignal';
% R = CSM(inputSignal, f0, fs, 64, 10);
R = R + trace(R)/(P^2)*eye(P, P);
R = R/P;
R = inv(R);

%Calculate weighting vector 
w = zeros(M, N, P);
for y = 1:M
    for x = 1:N
        ee = reshape(e(y, x, :), P, 1);
        wmvdr = R * ee / (ee' * R * ee);
        w(y, x, :) = wmvdr;
    end
end
% w = squeeze(w);
