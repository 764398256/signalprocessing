clear all;
close all;

% Config
power = 0.0001;

micInFile = '/home/mobvoi/work/data/kaiguanjiPC_ch3.wav';
spkInFile = '/home/mobvoi/work/data/kaiguanjiPC_ch4.wav';

[mic, mic_fs] = audioread(micInFile);
mic = mic(:,1);

[spk, spk_fs] = audioread(spkInFile);
spk = spk(:,1);

assert(mic_fs == spk_fs);
fs = mic_fs;

nSamples = min(size(spk, 1), size(mic,1));

frameSize = ceil(10 * fs / 1000);

nFrames = ceil(nSamples/frameSize);
spkPower = zeros(nFrames, 1);

for ii = 1:nFrames-1
   spkPower(ii) = var(spk((ii-1)*frameSize + (1:frameSize)));
end;

vad = zeros(nSamples, 1);

factor = 0.8;
for ii = 1:frameSize:nSamples
    if spkPower(ceil(ii/frameSize)) > power
        temp = 1;
    else
        temp = 0;
    end;
    if ii == 1
        vad(ii + (0:frameSize-1)) = temp;
    else
        vad(ii + (0:frameSize-1)) = factor * vad(ii -1) + (1 - factor) * temp;
    end;
end;

vad(vad > 0.01) = 1;
vad(vad ~= 1) = 0;

mic = mic.*mic;
spk = spk.*spk;    
mic = smooth(mic, 2);
spk = smooth(spk, 2);

riseIndex = 0;
descentIndex = 0;

steps = zeros(1, 1);
jj = 1;
for ii = 2:nSamples
    if(vad(ii-1) > 0 && vad(ii) == 0)
        descentIndex = ii;
    end;
    if(vad(ii-1) == 0 && vad(ii) > 0 && descentIndex > 0)
        riseIndex = ii;
        steps(jj) = ceil((descentIndex + riseIndex)/2);
        jj = jj + 1;
        riseIndex = 0;
        descentIndex = 0;
    end;
end;

assert(length(steps) > 1);

preDelay = 0;

delayTime = zeros(nSamples, 1);
for ii = 2:length(steps)
    len = steps(ii) - steps(ii-1);
    if len < fs
        delayTime(steps(ii-1)) = preDelay;
    else
        A = xcorr(spk(steps(ii-1):steps(ii), 1), mic(steps(ii-1):steps(ii), 1));
        [value, index] = max(A);
        delay = (steps(ii) -  steps(ii-1)) - index;
        delayTime(steps(ii-1)) = delay / fs;
        preDelay = delay / fs;
    end;
end;

figure(100);
time = (0:length(delayTime)-1)/fs;
plot(time, delayTime);