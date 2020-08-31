function w = weightingVectorLCMV(inputSignal)
%weighting Vector Linearly Constrained Minimum Variance 
%
%Calculates the Linearly Constrained Minimum Variance  for a specific array with
%a specific input signal
%
%w = weightingVectorLSMV(inputSignal)
%
%IN
%inputSignal         - PxL vector of inputsignals consisting of L samples

%
%OUT
%w                   - 1xP vector of complex optimal weights
%

%M # of theta steering angles, N # of phi steering angless, P number of mics
[P, ~] = size(inputSignal);

%Calculate correlation matrix with diagonal loading
R = inputSignal*inputSignal';
% R = CSM(inputSignal, f0, fs, 64, 10);

%Derive norm vector
u = zeros(P, 1);
if ~exist('normidx', 'var')
    normidx = ceil(P/2);
end
u(normidx) = 1;
uH = u';

%Calculate weighting vector 
invR = R/uH;
w = invR / (uH * invR) ;
