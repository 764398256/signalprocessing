clear all;
close all;

% Config
dtmf = 0;

micInFile = '/home/mobvoi/work/data/13_mic.wav';
spkInFile = '/home/mobvoi/work/data/13_ref.wav';

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

mic = smooth(mic, 2);

spk = smooth(spk, 2);

delayTime = zeros(nSamples, 1);

if dtmf == 1    
    nFrames = nSamples;
    vad = zeros(nFrames, 1);
    factor = 0.99;
    power = 0.01;
    for ii = 1:nFrames
        if spk(ii, 1)^2 > power
            temp = 1;
        else
            temp = 0;
        end;
        if ii == 1
            vad(ii) = temp;
        else
            vad(ii) = factor * vad(ii -1) + (1 - factor) * temp;
        end;
    end;
    
    if(0)
        figure(99);
        plot(spk, 'g');
        hold on;
        plot(vad* 0.1, 'r');
    end;
    
    startIndex = 0;
    stopIndex = 0;
    
    offset = 1600;
    
    for ii = 1:nFrames
        if(vad(ii) > 0.001 && startIndex == 0)
            startIndex = ii;
        elseif(vad(ii) < 0.001 && startIndex > 0)
            stopIndex = -1;
        elseif(vad(ii) > 0.001 && stopIndex == -1)
            stopIndex = ii;
            startIndex = startIndex - offset;
            stopIndex = stopIndex - offset;
            A = xcorr(spk(startIndex:stopIndex, 1), mic(startIndex:stopIndex, 1));
            [value, index] = max(A);
            delay = (stopIndex - startIndex) - index;
            delayTime(ii) = delay / 16000;
            startIndex = 0;
            stopIndex = 0;
        end;
    end;
    
else
    
    stepLen = ceil(500 * fs / 1000);
    overlapLen = ceil(stepLen / 2);
    delayTime = zeros(nSamples, 1);
    for ii = 1:overlapLen:nSamples-stepLen
        A = spk(ii + (1:stepLen));
        B = mic(ii + (1:stepLen));
        result = xcorr(A, B);
        [value, index] = max(result);
        delay = stepLen - index;
        delayTime(ii) = delay / fs;
    end;
end;

figure(100);
plot(delayTime);