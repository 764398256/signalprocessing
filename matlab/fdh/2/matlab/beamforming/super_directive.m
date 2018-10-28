clear all

f = 1000;
c = 340;
d = c/(2*f);  %0.17
% d = 0.05;
P = 9;
MArray = [-4*d,-3*d,-2*d,-d,0,d,2*d,3*d,4*d]';

%% 传统beamforming导向矢量
theta = 0;
A = exp(-j*2*pi*f*MArray*sind(theta)/c);
win = ones(length(A),1)/length(A);
w_beam = A.*win;

%% 超指向beamforming导向矢量
gama = zeros(P,P);
 for ch_i = 1:P
        for ch_j = ch_i:P
            dis = norm(MArray(ch_i)-MArray(ch_j));
            gama(ch_i,ch_j) = sinc(2*f*dis/c);
            gama(ch_j,ch_i) = gama(ch_i,ch_j);
        end
 end
As = exp(-j*2*pi*f*MArray*sind(theta)/c);
w_super= inv(gama)*As/(As'*inv(gama)*As);
   
%% 
J = zeros(181,1);
J_super = J;
idx = 1;
for theta = -90:1:90
    a = exp(-j*2*pi*f*MArray*sind(theta)/c);
    J(idx) = abs(w_beam'*a);
    J_super(idx) = abs(w_super'*a);
    idx = idx+1;
end

J = 20*log(J);
J_super = 20*log(J_super);
plot(-90:1:90,J)
hold on
plot(-90:1:90,J_super,'r');

