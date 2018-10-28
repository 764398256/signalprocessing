clear all
fs = 16000;
x = wavread('speech_noise.wav');
[B, f, t] = spectrogram(x./32768,512,256,512,fs);
imagesc(t,f,20*log10(abs(B))); 
title('speech+noise');
axis xy;

cohen = wavread('speech_cohen.wav');
figure
[B, f, t] = spectrogram(cohen./32768,512,256,512,fs);
imagesc(t,f,20*log10(abs(B))); 
title('Cohen');
axis xy;

webrtc = pcmread('speech_webrtc.pcm');
figure
[B, f, t] = spectrogram(webrtc./32768,512,256,512,fs);
imagesc(t,f,20*log10(abs(B))); 
title('Webrtc');
axis xy;

speex = wavread('speech_speex.wav');
figure
[B, f, t] = spectrogram(speex./32768,512,256,512,fs);
imagesc(t,f,20*log10(abs(B))); 
title('Speex');
axis xy;