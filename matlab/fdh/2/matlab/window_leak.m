clear all

f = 1000;
fs = 5120;
n = 1:512;
x = sin(2*pi*f*n/fs)+1*sin(2*pi*(f-10)*n/fs);
plot(x);
axis([1,512,-2,2])
X = abs(fft(x))/256;

xx = 0:10:2560;
plot(xx,X(1:257))
xlabel('ÆµÂÊ');
ylabel('·ù¶È');