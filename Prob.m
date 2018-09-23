close all;clear all;clc;

inFile = '/home/mobvoi/work/github/dsp/data/combined/f2_all.wav';
[mic, fs] = audioread(inFile);

[A, voiceFlags] = SnrVad(mic, fs);
mic = mic(voiceFlags > 0);

frameSize = 100 * fs / 1000;
Nh = 2 * frameSize;
% buffer x into a 2D matrix
x2 = buffer(mic, Nh, Nh - frameSize);
[N1, N2] = size(x2);
win = hanning(Nh);
% applying window funciton
x2 = x2.*repmat(win, 1, N2);
mic = reshape(x2, N1*N2, 1);

N = length(mic);
mic = mic/max(abs(mic));
M = 100;
prob = zeros(M, 1);
segm = linspace(-1,1,M+1);
for ii = 1:M
    prob(ii) = length(find(mic >= segm(ii) & mic < segm(ii+1)))/N;
end;
expect = mean(mic);
variance = var(mic);
norm = normpdf(segm(1:M), expect, sqrt(variance));

figure(1);
plot(segm(1:M), db(prob)/2, 'r');hold on;
plot(segm(1:M), db(norm/50)/2, 'g');