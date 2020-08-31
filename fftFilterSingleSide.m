function y=fftFilterSingleSide(y,Fs,fl,fu,guard)
flag = 0;
if size(y,2) ~= 1
    y = y.';
    flag = 1;
end
ny=numel(y);
if (ny<Fs)
    y=[y;zeros(Fs-ny,1)];
end
Y=fft(y);
N=size(Y,1);
Nl=floor(fl/Fs*N);
Nu=ceil(fu/Fs*N);
Ng=floor(guard/Fs*N);

filter=zeros(size(Y,1), size(Y,2));
filter(Nl:Nu)=1;
filter(Nl-Ng:Nl)=cos((Ng:-1:0)/Ng*pi/2);
filter(Nu:Nu+Ng)=cos((0:Ng)/Ng*pi/2);
% filter(N+2-Nu:N+2-Nl)=1;
% filter(N+2-Nl:N+2-Nl+Ng)=cos((0:Ng)/Ng*pi/2);
% filter(N+2-Nu-Ng:N+2-Nu)=cos((Ng:-1:0)/Ng*pi/2);
Y=Y.*filter;
y=ifft(Y);
y=y(1:ny);
if flag
    y = y.';
end