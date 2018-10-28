clear all
Fs=5120;
n=0:1/Fs:1;
f0=1000;
xn=sin(2*pi*f0*n)+sin(2*pi*(f0+19)*n)+0.1*randn(size(n));
window=boxcar(length(xn));
nfft=512;
[Pxx,f]=periodogram(xn,window,nfft,Fs);
figure(1)
plot(f,10*log10(Pxx)),grid
axis([0,2560,min(10*log10(Pxx)),max(10*log10(Pxx))])
xlabel('Frequency(Hz)')
ylabel('Power Spectral Density(dB/Hz)')
title('Periodogram PSD Estimate')
order1=20;
order2=50;
figure(2)
pburg(xn,order1,nfft,Fs)
figure(3)
pburg(xn,order2,nfft,Fs)