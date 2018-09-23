%clear all;close all;clc;

% Config
deltaTime = 0.6;
delta = 6;

melLen = 26;

testcase = 1;
if testcase == 1
    inFile = '/home/mobvoi/work/github/dsp/data/combined/f2_all.wav';
end;

[mic, fs] = audioread(inFile);
mic = mic(:,1);
assert(fs == 16000);

[A, voiceFlags] = SnrVad(mic, fs);
mic = mic(voiceFlags > 0);

frameSize = 10 * fs / 1000;
Nh = 3 * frameSize;
% buffer x into a 2D matrix
x2 = buffer(mic, Nh, Nh - frameSize);
[N1, N2] = size(x2);
win = hanning(Nh);
% applying window funciton
x2 = x2.*repmat(win, 1, N2);

MIC = fft(x2);
MIC = abs(MIC);

nfft = N1;

assert(fs == 16000);
startFreq = 100;
endFreq = 7800;
startMel= 1125*log(1+startFreq/700);
endMel = 1125*log(1+endFreq/700);
mels = linspace(startMel,endMel, melLen + 2);
freqs = 700*(exp(mels/1125)-1);
freqBins = floor(freqs*nfft/fs)+1;

coeff = zeros(melLen, N2);


for ii = 2:melLen +1
    temp = zeros(N1, 1);
    for iFreq = 1:ceil(freqBins(end)*fs/nfft);
        if(iFreq > freqBins(ii-1) && iFreq <= freqBins(ii))
            temp(iFreq) = (iFreq - freqBins(ii-1))/(freqBins(ii) - freqBins(ii-1));
        elseif(iFreq > freqBins(ii) && iFreq <= freqBins(ii+1))
            temp(iFreq) = (-iFreq + freqBins(ii+1))/(freqBins(ii+1) - freqBins(ii));
        end;
    end
    coeff(ii-1, :) = sum(MIC.*repmat(temp, 1,N2), 1);
end;

coeff = log(coeff);
coeff = dct(coeff);

figure(103);
for ii = 1:N2
    plot(coeff(:,ii)); hold on;
end
