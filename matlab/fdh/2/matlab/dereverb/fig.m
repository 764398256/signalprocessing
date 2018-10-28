clear all

[x,fs] = wavread('reverb_signal.wav');
[B, f, t] = spectrogram(x,512,256,512,fs);
imagesc(t,f,20*log10(abs(B))); 
axis xy;

wpe = wavread('dereverb_wpe');
figure
[B, f, t] = spectrogram(wpe,512,256,512,fs);
imagesc(t,f,20*log10(abs(B))); 
axis xy;

hab = wavread('dereverb_habits.wav');
figure
[B, f, t] = spectrogram(hab,512,256,512,fs);
imagesc(t,f,20*log10(abs(B))); 
axis xy;