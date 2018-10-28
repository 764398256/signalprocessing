clear all

x = wavread('SA1.wav');
fs = 16000;
t = 1/fs:1/fs:length(x)/fs;
plot(t,x);
[B, f, t] = spectrogram(x,256,128,256,fs);
figure
imagesc(t,f,20*log10(abs(B))); 
axis xy;

x = x(0.85*fs:fs);

win1 = 0*ones(length(x),1);
win2 = win1;
win1(1:512) = win;
win2(256:256+512-1) = win;
plot(x*3);
hold on
plot(win1,'r');
plot(win2,'r');
axis([1,length(x),-0.6,1]);
figure
plot(win1+win2,'r')
axis([1,length(x),-0.6,1.1]);