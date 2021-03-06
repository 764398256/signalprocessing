clear all

global  P   d  c M  f0

P = 8;
c = 340;
M = 2;
f0 = 500;
d = c/(2*f0);%circle array radius��
d = 0.17
snapshot = 100;
SNR = 10;
theta0 = [-10,5]*pi/180;
%% ���ɷ����ź�
coefficient = j*2*pi;
x = zeros(P,snapshot);

for i = 1:P
    for m= 1:M
        tao(i,m) = (i-1)*d*sin(theta0(m))/c;
    end
end
A = exp(coefficient*f0*tao);

s = randn(M,snapshot)+j*randn(M,snapshot);
sigma = 1/(10^(SNR/10));
w = sigma*(randn(P,snapshot)+randn(P,snapshot));
x = A*s+w;
%% 
Rx = zeros(P,P);
for n = 1:snapshot
    Rx = Rx+x(:,n)*x(:,n)';
end
Rx = Rx./snapshot;

step = 0.1;
Jc = zeros(1,180/step);

t = 0;
for k = 1:step:180.9
    t = t+1;
    theta = (k-90);
    Jc(t) = conventional_DOA(Rx,theta);
end

y = zeros(1,180/step);
theta0 = theta0*180/pi+89*ones(1,M);
theta0 = theta0/step+ones(1,M);
for i = 1:M
    y(theta0(i)) = 1;
end
xx = -89:step:90.9;

hold on
plot(xx,10*log10(Jc./max(max(Jc))));

plot(xx,10*log10(y+1e-10));

% axis([-60,60,-40,1]);
axis([-90,90,-40,1]);
xlabel('Degree')
ylabel('dB')



