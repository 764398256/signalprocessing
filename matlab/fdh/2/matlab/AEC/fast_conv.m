clear all

load('wo.mat');
plot(wo);

L = 512;
h = wo(1:L);
s = randn(L*20,1);

M = 128;
x = conv(s,h);
plot(x);

%% method 1
W = zeros(2*M,L/M);
for i = 1:L/M
    tmp = h((i-1)*M+1:i*M);
    W(1:M,i) = tmp;
    W(:,i) = fft(W(:,i));
end
n_blocks = length(s)/M;
Y = zeros(2*M,L/M);
in = zeros(2*M,1);
y = zeros(length(s),1);
for frame_n = 1:n_blocks
    in = zeros(2*M,1);
    in(1:M) = s((frame_n-1)*M+1:frame_n*M);
    S = fft(in);
    Y = [S,Y(:,1:end-1)];
    tmp = sum(Y.*W,2);
    tmp = ifft(tmp);
    if frame_n == 1
        y((frame_n-1)*M+1:frame_n*M) = tmp(1:M);
    else
        y((frame_n-1)*M+1:frame_n*M) = old+tmp(1:M);
    end
    old = tmp(M+1:end);
end

%%  method 2 
% W = zeros(2*M,L/M);
% for i = 1:L/M
%     tmp = h((i-1)*M+1:i*M);
%     W(1:M,i) = tmp;
%     W(:,i) = fft(W(:,i));
% end
% 
% n_blocks = length(s)/M;
% Y = zeros(2*M,L/M);
% in = zeros(2*M,1);
% y = zeros(length(s),1);
% for frame_n = 1:n_blocks
%     in(M+1:end) = s((frame_n-1)*M+1:frame_n*M);
%     S = fft(in);
%     Y = [S,Y(:,1:end-1)];
%     tmp = sum(Y.*W,2);
%     tmp = ifft(tmp);
%     y((frame_n-1)*M+1:frame_n*M) = tmp(M+1:end);
%     
%     in(1:M) = s((frame_n-1)*M+1:frame_n*M);
% end
%%

for i = 1:length(s)
    e(i) = sum(abs(x(i)-y(i)));
end
plot(e)
    
        

