clear alll

f = 1000;
c = 340;
d = c/(2*f);
P = 9;
MArray = [-4*d,-3*d,-2*d,-d,0,d,2*d,3*d,4*d]';

theta_s = 0;
As = exp(-j*2*pi*f*MArray*sind(theta_s)/c);

theta_i = [-50,30];
Ai = exp(-j*2*pi*f*MArray*sind(theta_i)/c);

Snap = 100;
ss = randn(1,Snap);

sigma_i = 20;
si = sigma_i*randn(2,Snap);
sigma_n = 1;
x = As*ss+Ai*si+sigma_n*randn(P,Snap);

Rx = zeros(P,P);
for i = 1:P
    Rx = Rx+x(:,i)*x(:,i)';
end
Rx = Rx./Snap;

sigma_l = 10;
Rx = Rx+sigma_l*eye(P);

A = inv(Rx)*As/(As'*inv(Rx)*As);
J = zeros(181,1);
idx = 1;
for theta = -90:1:90
    a = exp(-j*2*pi*f*MArray*sind(theta)/c);
    J(idx) = abs(a'*A);
    idx = idx+1;
end

% J = J./max(J);
J = 20*log(J);
plot(-90:1:90,J)
% axis([-90,90,-200,0])