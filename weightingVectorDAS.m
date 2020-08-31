function w = weightingVectorDAS(e)
%weightingVectorDAS - Delay and Sum Beamformer
%
%Calculates the Delay and Sum Beamformer weights for a specific array with
%a specific input signal and specific scanning angle
%
%w = weightingVectorMinimumVariance(inputSignal, e)
%
%IN
%e                   - 1xP steering vector of P antennas

%
%OUT
%w                   - 1xP vector of complex optimal weights
%

%M # of theta steering angles, N # of phi steering angless, P number of mics
[M, N, P] = size(e);

%Calculate weighting vector 
w = zeros(M, N, P);
for y = 1:M
    for x = 1:N
        ee = reshape(e(y, x, :), P, 1);
        w(y, x, :) = ee;
    end
end
% w = squeeze(w);

end