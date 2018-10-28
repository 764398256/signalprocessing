function [V_out,D_out] = RSS_eig(X,theta,f0)
% Circle array music J function
% theta
% N              the length of data
% R              circle array senser number
% Fs             signal frequency
% r              circle array's radius
% c              velocity of sound
% L              the number of samples in a segment
global  P   d  c  Fl Fh

bin_num = size(X,3);
deltf = (Fh-Fl)/(bin_num-1);
coefficient = -j*2*pi;

tao = zeros(P,length(theta));
for i = 1:P
    for m = 1:length(theta)
        tao(i,m) =(i-1)*d*sin(theta(m))/c;
    end
end

K = size(X,2);
RX = zeros(P,P,bin_num);
Rt = zeros(P,P);
A0 = exp(coefficient*f0*tao);
for n = 1:bin_num
    f = Fl+deltf*(n-1);
    for k = 1:K
        RX(:,:,n) = RX(:,:,n)+X(:,k,n)*X(:,k,n)';
    end
    RX(:,:,n) = RX(:,:,n)./K;
    A = exp(coefficient*f*tao);
    B = A*A0';
    [U,S,V] = svd(B);
    T = V*U';
    Rt = Rt+T*RX(:,:,n)*T';
end
[V_out,D_out] = eig_sort(Rt);
