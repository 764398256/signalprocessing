clear all

[data,fs] = wavread('SA1.wav');
output = zeros(length(data),1);
%% hanning window overlap
NFFT = 512;
overlap = 0.5;
win = hanning(NFFT);
frmshift = NFFT*overlap;
n_blocks = floor(length(data)/frmshift);
xbuf = zeros(NFFT,1);
outBuf = zeros(frmshift,1);
for n = 1:n_blocks
    xbuf(frmshift+1:end) = data((n-1)*frmshift+1:n*frmshift);
    X = fft(xbuf.*win);
    Y= [];
    for FreIdx = 1:NFFT/2+1
        Y(FreIdx,1) = X(FreIdx);
        if mod(FreIdx,3) == 0
            Y(FreIdx,1) = X(FreIdx)*2;
        end
    end
    
    Y = [Y; flipud(conj(Y(2:end-1)))];
    yinp = real(ifft(Y));  
   
    output((n-1)*frmshift+1:n*frmshift) = outBuf +  yinp(1 : frmshift);
    outBuf = yinp(frmshift+1:end);
    
    xbuf(1:frmshift,:) = xbuf(frmshift+1:end,:);
end;
wavwrite(output,fs,'SA1_process_overlap.wav');

%% no overlap
% NFFT = 512;
% win = hanning(NFFT);
% win = ones(NFFT,1);
% n_blocks = floor(length(data)/NFFT);
% xbuf = zeros(NFFT,1);
% frmshift = NFFT;
% for n = 1:n_blocks
%     xbuf = data((n-1)*frmshift+1:n*frmshift);
%     X = fft(xbuf.*win);
%     Y= [];
%     for FreIdx = 1:NFFT/2+1
%         Y(FreIdx,1) = X(FreIdx);
%         if mod(FreIdx,3) == 0
%             Y(FreIdx,1) = X(FreIdx)*2;
%         end
%     end
%     
%     Y = [Y; flipud(conj(Y(2:end-1)))];
%     yinp = real(ifft(Y));  
%    
%     output((n-1)*frmshift+1:n*frmshift) = yinp;
% end;
% wavwrite(output,fs,'SA1_process.wav');