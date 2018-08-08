clear all

% STFT parameters
windowLen = 400;
overlapLen = 160;
fftLen = 512;

% AEC parameters
echoLen  = 160;     % in ms
frameLen  = 10;     % in ms
nu_nlms = 0.8;

% Configs
algorithm = 2;
testcase = 2;

if testcase == 1
    micInFile = 'data/noise_mic20.wav';
    spkInFile = 'data/noise_spk.wav';
    micOutFile = 'data/out.wav';
elseif testcase == 2
    micInFile = 'data/micin_16k.wav';
    spkInFile = 'data/spkin_16k.wav';
    micOutFile = 'data/out.wav';
elseif testcase == 3
    micInFile = 'data/mic_c0.wav';
    spkInFile = 'data/spk.wav';
    micOutFile = 'data/out.wav';
end;
assert(frameLen == 10 || frameLen == 20);

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

if (fs == 16000)
    freqStart = 100;
    freqEnd = 7800;
    lowFreq = 1000;
elseif (fs == 48000)
    freqStart = 100;
    freqEnd = 14000;
    lowFreq = 1000;
else
    error('unsupported sampling rate');
end;

assert(fs/2 > freqEnd);

frameSize = frameLen*fs/1000;
numTaps = ceil(echoLen/frameLen);
numBands = frameSize * 2;
bandInterval = 1000/frameLen/2;    % 2x over-sampling filterbank

firstBand = freqStart/bandInterval;       % C style index
lastBand  = freqEnd/bandInterval - 1;     % C style index
usedBands = (firstBand+1:lastBand+1);

[micSpec, f, t] = Stft(mic, windowLen, overlapLen, fftLen, fs);
[spkSpec, f, t] = Stft(spk, windowLen, overlapLen, fftLen, fs);

if algorithm == 0
    [cleSpec, erle] = NlmsLi(micSpec, spkSpec, numTaps, nu_nlms);
elseif algorithm == 1
    cleSpec = Nlms(micSpec, spkSpec, numTaps, nu_nlms);
elseif algorithm == 2
    cleSpec = Apa(micSpec, spkSpec, numTaps, nu_nlms);
end;

[cleOut, t] = Istft(cleSpec, windowLen, overlapLen, fftLen, fs);

% Inverse filterbank analysis

audiowrite(micOutFile, cleOut, fs);