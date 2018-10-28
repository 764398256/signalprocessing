function [V_out,D_out] = TCT_f0_eig(X,theta,f0)

global  P   d  c M Fl Fh

coefficient = -j*2*pi;

tao = zeros(P,length(theta));
for i = 1:P
    for m = 1:length(theta)
        tao(i,m) =(i-1)*d*sin(theta(m))/c;
    end
end

K = size(X,2); 
bin_num = size(X,3);
deltf = (Fh-Fl)/(bin_num-1);
RX = zeros(P,P,bin_num);
V = RX;
S = zeros(length(theta),length(theta));
count = 0;
for n = 1:bin_num
    count = count+1;
    f = Fl+deltf*(n-1);
    A = exp(coefficient*f*tao);
    for k = 1:K
        RX(:,:,n) = RX(:,:,n)+X(:,k,n)*X(:,k,n)';
    end
    RX(:,:,n) = RX(:,:,n)./K;
    [v,dj] = eig_sort(RX(:,:,n));
    Pj = RX(:,:,n)-min(dj)*eye(P);
    [V(:,:,n),D(n,:)] = eig_sort(Pj);
    B = inv(A'*A)*A';
    S = S+B*Pj*B';
end
S = S/(count);
%Çóf0
% mu = zeros(1,M);
% for i = 1:M
%     count = 0;
%     for n = 1:bin_num
%         count = count+1;
%         mu(i) = mu(i)+D(n,i);
%     end
% end
% mu = mu/count;
% f_min = 1e+30*ones(1,bin_num);
% for n = 1:bin_num
%     f = Fl+deltf*(n-1);
%     A0 = exp(coefficient*f*tao);
%     P0 = A0*S*A0';
%     [V0,D0] = eig_sort(P0);
%     f_min(n) = 0;
%     for i = 1:M
%         f_min(n) = f_min(n)+(abs(D0(i)-mu(i)))^2;
%     end
% end
% f0 = find(f_min == min(min(f_min)));
% f0 = Fl+deltf*(f0-1);
A0 = exp(coefficient*f0*tao);
P0 = A0*S*A0';
[V0,D0] = eig_sort(P0);
RX_a = zeros(P,P);
count = 0;
for n = 1:bin_num
    count = count+1;
    U = V0*V(:,:,n)';
    RX_a = RX_a+U*RX(:,:,n)*U';
end
RX_a = RX_a./(count);
[V_out,D_out] = eig_sort(RX_a);
