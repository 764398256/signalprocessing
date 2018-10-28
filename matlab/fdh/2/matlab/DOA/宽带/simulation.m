clear all;% 

global  P   d  c M Fl Fh

P = 8;
c = 340;

M = 2;
Fl = 80;
Fh = 120;
f0 = 100;
bin_num = 9;
d = c/(2*f0);%circle array radiusª¤
snapshot = 50;
SNR = 5;

amplitude = [1,1];
theta0 = [11,13]*pi/180;

tic;%¼ÆÊ±
X = generate_signal_Wbins(Fl,Fh,bin_num,snapshot,SNR,amplitude,theta0);

% [V,D] = music_eig(X);

% BW = 0.886*lamda/(R*d) 6.4
theta = [9,12,15]*pi/180;
% theta = [10,11.5,14]*pi/180;
% theta = [11,13]*pi/180;

[V_r,D_r] = RSS_eig(X,theta,f0);

[V_T,D_T,f0_T] = TCT_eig(X,theta);
% [V_T,D_T] = TCT_f0_eig(X,theta,f0);

step = 0.1;
Jm = zeros(1,180/step);
Jr = Jm;
Jt = Jm;
t = 0;
for k = 1:step:180.9
    t = t+1;
    theta = (k-90)*pi/180;
%     Jm(t) = music_DOA(V,theta);
    Jr(t) = RSS_DOA(V_r,f0,theta);
%     Jt(t) = TCT_DOA(V_T,f0,theta);
    Jt(t) = TCT_DOA(V_T,f0_T,theta);

end

y = zeros(1,180/step);
theta0 = theta0*180/pi+89*ones(1,M);
theta0 = theta0/step+ones(1,M);
for i = 1:M
    y(theta0(i)) = 1;
end

xx = -89:0.1:90.9;
% plot(xx,10*log10(Jm./max(max(Jm))))
hold on

plot(xx,10*log10(Jr./max(max(Jr))),'g')

plot(xx,10*log10(Jt./max(max(Jt))),'r')

plot(xx,10*log10(y+1e-10),'--');

axis([-20,40,-40,1]);
xlabel('Degree')
ylabel('dB')
% legend('ISM','ISM-Smooth','RSS','RSS-Smooth','TCT','TOPS')



ave = toc



