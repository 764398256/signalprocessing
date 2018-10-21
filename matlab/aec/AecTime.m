clear all;clc;close all;

testcase = 2;

if testcase == 1
    micInFile = 'data/noise_mic20.wav';
    spkInFile = 'data/noise_mic20.wav';
    micOutFile = 'data/out.wav';
elseif testcase == 2
    micInFile = 'data/micin_16k.wav';
    spkInFile = 'data/micin_16k.wav';
    micOutFile = 'data/out.wav';
end;

[mic, mic_fs] = audioread(micInFile);
if size(mic, 2) > 1
    mic(:, 2:end) = [];
end;

[spk, spk_fs] = audioread(spkInFile);
if size(spk, 2) > 1
    spk(:, 2:end) = [];
end;

assert(mic_fs == spk_fs);      % does not support re-sampling
fs = mic_fs;

r = Rir(mic_fs, 0.1);
mic = conv(mic, r);

% Configs
nTaps = length(r);
nu_nlms = 0.8;

nSamples = min(size(spk, 1), size(mic, 1));
h = zeros(nTaps, 1);
clean = zeros(nSamples, 1);
for ii = nTaps:nSamples
    x = spk(ii-nTaps+1:ii);
    e = mic(ii, 1) - h.'*x;
    power = x.'*x;
    % if (power/nTaps > 0.032^2)
        h = h + nu_nlms*(e*x./power);
    % end;
    clean(ii) = e;
end;

frameSize = fs*10/1000;
windowSize = frameSize*2;

input2 = buffer(mic, windowSize, windowSize-frameSize, 'nodelay');
output2 = buffer(clean, windowSize, windowSize-frameSize, 'nodelay');

N2 = min([size(input2,2) size(output2,2)]);
input2 = input2(:, 1:N2);
output2 = output2(:, 1:N2);

erleSorted = sort(10*log10(var(input2)./var(output2)));
e95p = erleSorted(round(length(erleSorted)*0.95));

fprintf('erle is %f\n', e95p);

figure(90);
plot(erleSorted);
grid on;

figure(91);
plot(clean, 'b');
grid on;

audiowrite(micOutFile, clean, fs);
