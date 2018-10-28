clear alll

f = 1000;
c = 340;
d = c/(2*f);
P = 9;
MArray = [-4*d,-3*d,-2*d,-d,0,d,2*d,3*d,4*d]';

theta = 0;
A = exp(-j*2*pi*f*MArray*sind(theta)/c);
win = ones(length(A),1)/length(A);
win = hanning(length(MArray))/length(A);
win = hamming(length(MArray))/length(A);
Aw = A.*win;
J = zeros(181,1);
idx = 1;
for theta = -90:1:90
    a = exp(-j*2*pi*f*MArray*sind(theta)/c);
    J(idx) = abs(a'*Aw);
    idx = idx+1;
end

J = J./max(J);
J = 20*log(J);
plot(-90:1:90,J)
axis([-90,90,-200,0])