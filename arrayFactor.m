function AF = arrayFactor(w, e)
%arrayFactor - Calculate array factor of 1D, 2D or 3D array
%
%This matlab function calculates the array factor of a 1D, 2D or 3D array based
%on the position of the elements/sensors and the weight associated with
%each sensor. 
%
%AF = arrayFactor(w,e)
%
%IN
%w      - 1xP vector of element weights
%e      - MxNxP vector of steering vector

%
%OUT
%AF     - MxNxP Calculated array factor
%

%M # of theta steering angles, N # of phi steering angless, P number of mics
[M, N, P] = size(e);

g = repmat(reshape(w, 1, 1, P), M, N);
AF = abs(sum(conj(g) .* e, 3));

%Normalising
AF = abs(AF./max(max(AF)));

%
%                 N
%AF(theta, phi) = sum [ g_n * exp{jk(u*x_n + v*y_n + w*z_n)} ]
%                n=1
%
