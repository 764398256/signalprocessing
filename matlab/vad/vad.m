function [vadIndex, speechPower] = vad(fname)

[mic, fs] = audioread(fname);
if size(mic, 2) > 1
    mic = mean(mic, 2);
end
N = length(mic);
t = (0:N-1)/fs;

% parameters
frameLen = 10; % 10ms

if fs == 16000
    startFreq = 100;
    stopFreq = 7800;
elseif fs == 48000
    startFreq = 100;
    stopFreq = 14000;
else
    error('unsupported sample rate');
end

frameSize = fs * frameLen / 1000;
numBands = frameSize * 2;
bandInterval = fs/numBands;
firstBandUsed = startFreq/bandInterval;     % C style index
numBandUsed = (stopFreq-startFreq)/bandInterval;
lastBandUsed = stopFreq/bandInterval - 1;       % Inclusive, C style index

h = loadFbCoef(numBands);
energyFactor = 1/sum(h.^2);

% Forward filterbank analysis
X = DftFB(mic, frameSize, numBands, h, length(h)/numBands);

% zero out low bands and high bands. Bandpass 100-14000 Hz
X(1:firstBandUsed, :) = 0;
usedSpec = X(firstBandUsed + (1:numBandUsed), :);
micPower = abs(usedSpec.^2);
micEnergy = sum(micPower, 1)/frameSize * energyFactor;
signalLevel = sqrt(micEnergy);
micNoisePower = NoiseEst(usedSpec, fs, bandInterval, firstBandUsed, energyFactor);
micNoiseEnergy = sum(micNoisePower, 1)/frameSize * energyFactor;
micNoiseLevel = sqrt(micNoiseEnergy);

% a simple VAD
snrThreshold = 4; % 12 dB
voiceBand = [200/bandInterval : 2500/bandInterval]; %200 ~ 2500 Hz, C style index
signalPowerVoiceBand = micPower(voiceBand - firstBandUsed + 1, :)/frameSize * energyFactor;
noisePowerVoiceBand = micNoisePower(voiceBand - firstBandUsed + 1, :)/frameSize * energyFactor;
voiceFlag = sum(signalPowerVoiceBand, 1) > sum(noisePowerVoiceBand, 1) * snrThreshold * snrThreshold;
vadIndex = find(voiceFlag>0);

speechPower = mean(micEnergy(vadIndex));

speechPower = db(speechPower)/2;

fprintf('speech power %f\n', speechPower);


t_frame = (0:size(X,2)-1)*frameLen/1000;
close all;
figure(10);
plot(t_frame, voiceFlag * 0.1, 'r'); hold on;
plot(t_frame, signalLevel, 'g');

if(0)
figure(11);
plot(t_frame, db(signalLevel));
hold on
plot(t_frame, db(micNoiseLevel),'r');
plot(t_frame(vadIndex), db(micNoiseLevel(vadIndex)), 'g.');
grid on
end;
