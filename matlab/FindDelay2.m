clear all;
close all;

% Config
threshold = 0.7;
stepSize = 300; % in ms

micInFile = '/home/mobvoi/work/data/0815/08151232/in_ch3.wav';
spkInFile = '/home/mobvoi/work/data/0815/08151232/in_ch4.wav';

[mic, mic_fs] = audioread(micInFile);
if size(mic, 2) > 1
    mic(:, 2:end) = [];
end;

[spk, spk_fs] = audioread(spkInFile);
if size(spk, 2) > 1
    spk(:, 2:end) = [];
end;

assert(mic_fs == spk_fs);
fs = mic_fs;

nSamples = min(size(spk, 1), size(mic,1));

for ii = 1:nSamples
    mic(ii,1) = mic(ii,1)^2;
end;

for ii = 1:nSamples
    spk(ii,1) = spk(ii,1)^2;
end;

mic = smooth(mic, 3);
spk = smooth(spk, 3);

stepLen = ceil(stepSize * fs / 1000);
overlapLen = ceil(stepLen / 2);
coeff = zeros(nSamples, 1);
delayTime = zeros(nSamples, 1);

for ii = 1:overlapLen:nSamples-stepLen
    A = spk(ii + (1:stepLen));
    B = mic(ii + (1:stepLen));
    result = xcorr(A, B)/(norm(A)*norm(B));
    [value, index] = max(result);
    delay = stepLen - index;
    coeff(ii) = value;
    if(value > threshold)
        delayTime(ii) = delay / fs;
    else
        if(ii ~= 1)
            delayTime(ii) = delayTime(ii - overlapLen);
        end;
    end;
end;

figure(100);
time = (0:length(delayTime)-1)/fs;
plot(time, delayTime);