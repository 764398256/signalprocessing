clear all
fs = 16000;
x = wavread('mono.wav');
[B, f, t] = spectrogram(x,512,256,512,fs);
imagesc(t,f,20*log10(abs(B))); 
axis xy;

x1 = wavread('Carioid_first.wav');
figure
[B, f, t] = spectrogram(x1,512,256,512,fs);
imagesc(t,f,20*log10(abs(B))); 
title('first order');
axis xy;

x2 = wavread('Carioid_second.wav');
figure
[B, f, t] = spectrogram(x2,512,256,512,fs);
imagesc(t,f,20*log10(abs(B))); 
title('second order');
axis xy;