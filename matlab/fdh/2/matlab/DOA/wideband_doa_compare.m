clear all;% 

global  P   d  c M Fl Fh

P = 8;
c = 340;

M = 2;
Fl = 80;
Fh = 120;
f0 = 100;
bin_num = 9;
d = c/(2*f0);%circle array radius
snapshot = 50;
SNR = 5;

theta0 = [11,13]*pi/180;

global  P   d  c M Fl Fh

P = 8;
c = 340;

M = 2;
Fl = 80;
Fh = 120;
f0 = 100;
bin_num = 9;
d = c/(2*f0);
snapshot = 50;
SNR = 5;

theta0 = [10,14]*pi/180;

%% generate wideband array signal
f = Fl:(Fh-Fl)/(bin_num-1):Fh;
coefficient = -j*2*pi;
X = zeros(P,snapshot,length(f));
for k = 1:length(f)
    S = zeros(M,snapshot);
    for i = 1:M
        S(i,:) = (randn(1,snapshot)+j*randn(1,snapshot));
    end
    for i = 1:P
        for m = 1:M
            tao(i,m) = (i-1)*d*sin(theta0(m))/c;
        end
    end
    sigma = 1/(10^(SNR/20));
    A = exp(coefficient*f(k)*tao);
    for n = 1:snapshot
        W = sigma*(randn(P,1)+j*randn(P,1));
        X(:,n,k) = A*S(:,n)+W;
    end
end

%% 计算MUSIC算法各个频带的特征值和特征向量
K = size(X,2); 
bin_num = size(X,3);
for n = 1:bin_num
    RX1 = zeros(P,P);
    for k = 1:K
        RX1 = RX1+X(:,k,n)*X(:,k,n)';
    end
    RX1 = RX1./K;
    [V(:,:,n),D(:,n)] = eig_sort(RX1);
end
%%

% BW = 0.886*lamda/(R*d) 6.4
theta = [9,12,15]*pi/180;
[V_r,D_r] = RSS_eig(X,theta,f0);

step = 0.1;
Jm = zeros(1,180/step);
Jr = Jm;
t = 0;
for k = 1:step:180.9
    t = t+1;
    theta = (k-90)*pi/180;
    Jm(t) = wideband_music_DOA(V,theta);
    Jr(t) = RSS_DOA(V_r,f0,theta);
end

y = zeros(1,180/step);
theta0 = theta0*180/pi+89*ones(1,M);
theta0 = theta0/step+ones(1,M);
for i = 1:M
    y(theta0(i)) = 1;
end

xx = -89:0.1:90.9;
plot(xx,10*log10(Jm./max(max(Jm))))
hold on

plot(xx,10*log10(Jr./max(max(Jr))),'r')

plot(xx,10*log10(y+1e-10),'--');

axis([-20,40,-40,1]);
xlabel('Degree')
ylabel('dB')



