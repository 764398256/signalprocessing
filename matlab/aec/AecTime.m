clear all;
close all;

nu_nlms = 0.8;

testcase = 2;

if testcase == 1
    micInFile = 'data/noise_mic20.wav';
    spkInFile = 'data/noise_spk.wav';
    micOutFile = 'data/out.wav';
elseif testcase == 2
    micInFile = 'data/micin_16k.wav';
    spkInFile = 'data/spkin_16k.wav';
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

% 100 ms data
nTaps = ceil(10 * 10 * fs / 1000);

nSamples = min(size(spk, 1), size(mic, 1));

h = zeros(nTaps, 1);
clean = zeros(nSamples, 1);

for ii = nTaps:nSamples
    x = spk(ii-nTaps+1:ii);
    e = mic(ii, 1) - h.'*x;
    power = x.'*x;
    if (power/nTaps > 0.032^2)
        h = h + nu_nlms*(e*x./power);
    end;
    clean(ii) = e;
end;

if(0)
figure(2);
plot(clean, 'b');
end;

audiowrite(micOutFile, clean, fs);
