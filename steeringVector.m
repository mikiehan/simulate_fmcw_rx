function e = steeringVector(xPos, yPos, zPos, f, c, thetaScanAngles, phiScanAngles)
%steeringVector - calculate steering vector of array
%
%Calculates the steering vector for different scanning angles
%Theta is the elevation and is the normal incidence angle
%Phi is the azimuth, and is the angle in the XY-plane
%
% e = steeringVector(xPos, yPos, zPos, f, c, thetaScanAngles, phiScanAngles)
%
%IN
%xPos            - 1xP vector of x-positions [m]
%yPos            - 1xP vector of y-positions [m]
%zPos            - 1xP vector of z-positions [m]
%f               - Wave frequency [Hz]
%c               - Speed of sound [m/s]
%thetaScanAngles - 1xM vector or MxN matrix of theta scanning angles [degrees]
%phiScanAngles   - 1xN vector or MxN matrix of of phi scanning angles [degrees]
%
%OUT
%e               - MxNxP matrix of steering vectors
%

if ~isvector(xPos)
    error('X-positions of array elements must be a 1xP vector where P is number of elements')
end

if ~isvector(yPos)
    error('Y-positions of array elements must be a 1xP vector where P is number of elements')
end

if ~isvector(zPos)
    error('Z-positions of array elements must be a 1xP vector where P is number of elements')
end

%Convert angles to radians
thetaScanAngles = thetaScanAngles*pi/180;
phiScanAngles = phiScanAngles*pi/180;
 
%Wavenumber
k = 2*pi*f/c;

%M # of y-points, N # of x-points, P number of mics
M = numel(thetaScanAngles);
N = numel(phiScanAngles);
P = numel(xPos);
ant_posi = [xPos; yPos; zPos]';

%exp(-1i*2*pi*freq*tau);

%incidentdir = [-cosd(elang).*cosd(azang);...
%    -cosd(elang).*sind(azang);...
%    -sind(elang)];
%tau = pos.'*incidentdir/c;

e = zeros(M, N, P);
for y = 1:M
        for x = 1:N
           angleT = [sin(phiScanAngles(x)) * cos(thetaScanAngles(y)); ...
                     sin(phiScanAngles(x)) * sin(thetaScanAngles(y)); ...
                     cos(phiScanAngles(x))]; 
           e(y, x, :) = exp(1j * k * (ant_posi * angleT));
        end
end

%         %steeringVec2(mic, angle) = cos(angle - 2*pi * (mic-1)/8);
%x2 = radius*(cos(angle) - cos(j*((360.0)/(8.0) - angle))
%ebi = cos(k*x);

