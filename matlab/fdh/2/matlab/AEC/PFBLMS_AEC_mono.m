clear all
fs = 16000;

%% simulation
% x = wavread('farend.wav');
% nearend = wavread('SA1.wav');
% wo = wavread('room_ir.wav');
% 
% farend = conv(x,wo);
% d = zeros(length(farend),1);
% d = d+farend;
% d(3*fs:3*fs+length(nearend)-1) = d(3*fs:3*fs+length(nearend)-1)+nearend;
% wavwrite(d,16000,'simu.wav');

%% real signal
data = wavread('real2.wav');
% data = wavread('sai_real.wav');
x = data(:,2)*20;
d = data(:,1);
%%
N=1024; M=128; P=N/M;
wF=zeros(2*M,P);
xF=zeros(2*M,P);
mu=0.5;
e=d;
epost = d;
zeta=zeros(size(d));

%%ÂË²¨Æ÷Îó²î
wF_old = wF;

outBuf = zeros(M,1);
xBuf = zeros(2*M,1);
dBuf = zeros(2*M,1);
eBuf = zeros(2*M,1);
Sd = zeros(M+1,1);
Se = zeros(M+1,1);
Sx = zeros(M+1,1);
Sde = zeros(M+1,1);
Sxd = zeros(M+1,1);
Sxe = zeros(M+1,1);
alpha = 0.9;
frame_idx = 1;
for n=M+1:M:length(x)-M
    xF=[fft(x(n-M:n+M-1)) xF(:,[1:end-1])];
    yhat=ifft(sum((wF.*xF).').'); 
    yhat=real(yhat(M+1:end));
    E=d(n:n+M-1)-yhat;
    MU=mu*(sum((abs(xF).^2)')'+0.1).^(-1);
    EF=fft([zeros(M,1); E]);
    for p = 1:P
        wF(:,p) = wF(:,p)+MU.*EF.*conj(xF(:,p));
    end
    %ÂË²¨Æ÷Îó²î
    WE(frame_idx) = sum(sum(abs((wF-wF_old))));
    wF_old = wF;
    frame_idx = frame_idx+1;
    %
    waux=real(ifft(wF)); wF=fft([waux(1:M,:); zeros(M,P)]);
    e(n:n+M-1)=d(n:n+M-1)-yhat;
 %% post filter
    eBuf(M+1:end) = e(n:n+M-1);
    efw = fft(eBuf.*hanning(M*2));
    efw = efw(1:M+1);
    
    dBuf(M+1:end) = d(n:n+M-1);
    dfw = fft(dBuf.*hanning(M*2));
    dfw = dfw(1:M+1);
    
    xBuf(M+1:end) = x(n:n+M-1);
    xfw = fft(xBuf.*hanning(M*2));
    xfw = xfw(1:M+1);
    
    Sd = alpha*Sd+(1-alpha)*(dfw.*conj(dfw));
    Se = alpha*Se+(1-alpha)*(efw.*conj(efw));
    Sx = alpha*Sx+(1-alpha)*(xfw.*conj(xfw));
    
    Sde = alpha*Sde+(1-alpha)*(dfw.*conj(efw));
    Sxe = alpha*Sxe+(1-alpha)*(xfw.*conj(efw));
    Sxd = alpha*Sxd+(1-alpha)*(xfw.*conj(dfw));
    
    cohde = Sde.*conj(Sde)./(Sd.*Se+1e-10);
    cohde = min(cohde,1);
    cohxd = Sxd.*conj(Sxd)./(Sx.*Sd+1e-10);
    cohxd = min(cohxd,1);
    
    hNl = min(cohde,1-cohxd);
    efw = efw.*hNl;
    ft = [efw;flipud(conj(efw(2:end-1,:)))];
    ft = ifft(ft);
    epost(n:n+M-1) = ft(1:M)+outBuf;
    outBuf = ft(M+1:end);
    
    xBuf(1:M) = xBuf(M+1:end);
    dBuf(1:M) = dBuf(M+1:end);
    eBuf(1:M) = eBuf(M+1:end);
end
% e = e(1:fs*10);
% epost = epost(1:fs*10);
% d = d(1:fs*10);
[B, f, t] = spectrogram(d,1024,512,512,fs);
figure
imagesc(t,f,20*log10(abs(B))); 
axis xy;
[B, f, t] = spectrogram(e,1024,512,512,fs);
figure
imagesc(t,f,20*log10(abs(B))); 
axis xy;
[B, f, t] = spectrogram(epost,1024,512,512,fs);
figure
imagesc(t,f,20*log10(abs(B))); 
axis xy;
wavwrite(e,fs,'aec_out.wav');
wavwrite(epost,fs,'aec_post_out.wav');

% xx = pcmread('wav/real/aec_out.pcm');
% [B, f, t] = spectrogram(xx,1024,512,512,fs);
% figure
% imagesc(t,f,20*log10(abs(B))); 
% axis xy;

