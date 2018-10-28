clear all
fs = 16000;

MArray = 0.035*[1,0;0,1;-1,0;0,-1];
MicNum = 4;
fs = 16000;
NFFT = 512;
c = 340;
data = wavread('sai_4mic.wav');

FreBand = linspace(0,fs/2,NFFT/2+1);
AnaWnd = hanning(NFFT);

FrmShift = NFFT/2;
n_blocks = floor(length(data)/FrmShift);
xbuf = zeros(NFFT,MicNum);
Xbuf = zeros(NFFT,MicNum);
DesAng = 0;

%% 传统波束形成加权向量
w_con = zeros(MicNum,NFFT/2+1);
for FreIdx = 1 : length(FreBand)
    coefficient = 2 * pi * FreBand(FreIdx)/c;
    h = [cosd(DesAng),sind(DesAng)].';
    tao = MArray*h;
    w_con(:,FreIdx) = exp(j*coefficient*tao)/MicNum;
end;
Output = zeros(length(data),1);
outBuf = zeros(NFFT/2,1);

%% 超指向性波束形成加权向量
w_super = zeros(MicNum,NFFT/2+1);
As = zeros(MicNum,NFFT/2+1);
for FreIdx = 1 : length(FreBand)
    coefficient = 2 * pi * FreBand(FreIdx)/ c;
    h = [cosd(DesAng),sind(DesAng)].';
    tao = MArray*h;
    As(:,FreIdx) = exp(j*coefficient*tao);
end

gama = zeros(MicNum,MicNum);
for FreIdx = 2:length(FreBand)
    for ch_i = 1:MicNum
        for ch_j = ch_i:MicNum
            dis = norm(MArray(ch_i,:)-MArray(ch_j,:));
            gama(ch_i,ch_j,FreIdx) = sinc(2*FreBand(FreIdx)*dis/c);
            gama(ch_j,ch_i,FreIdx) = gama(ch_i,ch_j,FreIdx);
        end
    end
    w_super(:,FreIdx) = inv(gama(:,:,FreIdx))*As(:,FreIdx)/(As(:,FreIdx)'*inv(gama(:,:,FreIdx))*As(:,FreIdx));
end
w_super(:,1:4) = w_con(:,1:4);
%% mvdr算法参数
Rx = zeros(MicNum,MicNum,NFFT/2+1);
w_mvdr = zeros(MicNum,NFFT/2+1);
w_mvdr(:,1) = w_con(:,1);
w_mvdr(:,NFFT/2+1) = w_con(:,NFFT/2+1);

%% 
w_beam = w_con;
for n = 1:n_blocks
    for ch = 1:MicNum
        xbuf(NFFT/2+1:end,ch) = data((n-1)*NFFT/2+1:n*NFFT/2,ch);
        Xbuf(:,ch) = fft(xbuf(:,ch).*AnaWnd);
    end
    
    %% mvdr 算法
    Xv = Xbuf(1:NFFT/2+1,:);
    for FreIdx = 1:NFFT/2+1
        Rx(:,:,FreIdx) = Rx(:,:,FreIdx)+Xv(FreIdx,:)'*Xv(FreIdx,:);
    end
    if mod(n,10) == 0
        Rx = Rx./10;
        sigma_l = 0.0001;
        for FreIdx = 2:NFFT/2
            Rx_l = Rx(:,:,FreIdx)+sigma_l*eye(MicNum);
            w_mvdr(:,FreIdx) = inv(Rx_l)*As(:,FreIdx)/(As(:,FreIdx)'*inv(Rx_l)*As(:,FreIdx));
        end
        Rx = zeros(MicNum,MicNum,NFFT/2+1);
        
        w_beam = w_mvdr;
    end
   %%  
    
    Y_beam = [];
    for FreIdx = 1 : NFFT/2+1
        X = Xbuf(FreIdx,:).';
        Y_beam(FreIdx,1) = w_beam(:,FreIdx)'*X;
    end;
    
    Y_beam = [Y_beam; flipud(conj(Y_beam(2:end-1)))];
    yinp = real(ifft(Y_beam));  
   
    Output((n-1)*NFFT/2+1:n*NFFT/2) = outBuf +  yinp(1 : NFFT/2);
    outBuf = yinp(NFFT/2+1:end);
    
    
    xbuf(1:NFFT/2,:) = xbuf(NFFT/2+1:end,:);
end;
wavwrite(Output,fs,'beam1.wav');
[B, f, t] = spectrogram(data(:,1),1024,512,1024,fs);
imagesc(t,f,20*log10(abs(B))); 
axis xy;

% Output = wavread('sai_output.wav');
% figure
% [B, f, t] = spectrogram(Output,1024,512,1024,fs);
% imagesc(t,f,20*log10(abs(B))); 
% axis xy;


