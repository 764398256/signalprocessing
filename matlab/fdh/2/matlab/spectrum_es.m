clear all

f = 1000;
fs = 5120;
% n = 1:512;
n=0:1/fs:1;
x = sin(2*pi*(f-15)*n)+1*sin(2*pi*(f-25)*n)+0*randn(1,length(n));
X = abs(fft(x))/256;

% xm1 = x(1:256);
% xm2 = x(129:256+129-1);
% xm3 = x(257:end);
% Xm = abs(fft(xm1))+abs(fft(xm2))+abs(fft(xm3));
% Xm = Xm/3/128;
% xx = 0:20:2560;
% plot(xx,Xm(1:129));
% xlabel('频率');
% ylabel('幅度');

figure
xx = 0:10:2560;
plot(xx,X(1:257))
xlabel('频率');
ylabel('幅度');

figure;
[Pxx,fx]=periodogram(x,[],256,fs);
plot(fx,Pxx)

[Xar,fx] = pburg(x,20,256,fs);
% [Xar,fx] = pmusic(x,30,512,fs);
figure
plot(fx,Xar)