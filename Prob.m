% close all;clear all;clc;

% inFile = '/home/mobvoi/work/github/dsp/data/noise/n5.wav'; 
% inFile = '/home/mobvoi/work/data/loudao1_ch1.wav';
inFile = '/home/mobvoi/work/github/dsp/data/combined/f2_all.wav';
outFile = '/home/mobvoi/work/github/dsp/data/combined/f2_all_out.wav';

[mic, fs] = audioread(inFile);

% % Add reverb for clean speech
% c = 340;                    % Sound velocity (m/s)
% fs = 16000;                 % Sample frequency (samples/s)
% r = [3 3 3];              % Receiver position [x y z] (m)
% s = [0 3 3];              % Source position [x y z] (m)
% L = [6 6 6];                % Room dimensions [x y z] (m)
% beta = 0.6;                 % Reverberation time (s)
% n = 2048;                   % Number of samples
% 
% h = rir_generator(c, fs, r, s, L, beta, n);
% 
% mic = filter(h, 1, mic);
% 
% audiowrite(outFile, mic, fs);

h = rir(fs, 1);

mic = filter(h, 1, mic);

% audiowrite(outFile, mic, fs);

voiceFlags = SnrVad(mic, fs);
mic = mic(voiceFlags > 0);

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
plot(segm(1:M), prob, 'r');hold on;
plot(segm(1:M), norm/50, 'g');


