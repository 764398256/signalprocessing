clear all
fs = 16000;
[h,fh] = wavread('perth_city_hall_balcony_ir_edit.wav');
h = resample(h,16000,fh);

data = wavread('singing.wav');
data = resample(data,fs,48000);
% data = wavread('song.wav');
data_g = conv(data,h(:,1));
data_g = data_g./max(abs(data_g))*0.5;
wavwrite(data_g,fs,'song_hall.wav');

[B, f, t] = spectrogram(data,1024,512,1024,fs);
imagesc(t,f,20*log10(abs(B))); 
axis xy;
figure
[B, f, t] = spectrogram(data_g,1024,512,1024,fs);
imagesc(t,f,20*log10(abs(B))); 
axis xy;
