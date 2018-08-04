clear all;
close all;

testcase = 1;

if testcase == 1
    micInFile = '/home/mobvoi/work/data/16_mic.wav';
    spkInFile = '/home/mobvoi/work/data/16_ref.wav';
end;

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

nFrames = min(size(spk, 1), size(mic,1));

for ii = 1:nFrames-2
    mic(ii,1) = (mic(ii,1)^2 + mic(ii+1,1)^2 + mic(ii+2, 1)^2) / 3;
end;

for ii = 1:nFrames-2
    spk(ii,1) = (spk(ii,1)^2 + spk(ii+1,1)^2 + spk(ii+2, 1)^2) / 3;
end;

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

% figure(100);
% plot(spk, 'g');
% hold on;
% plot(vad* 0.1, 'r');

startIndex = 0;
stopIndex = 0;

offset = 1600;

delayintime = zeros(nFrames, 1);

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
        
%         figure(102);
%         plot(spk(startIndex:stopIndex, 1)); hold on;
%         plot(mic(startIndex:stopIndex, 1));
        
        % fprintf('%f\n', delay / 16000);
        delayintime(ii) = delay / 16000;  
        startIndex = 0;
        stopIndex = 0;
    end;
end;

figure(103);
plot(delayintime);



