function J = MUSIC_DOA(V,theta)

global  P   d  c M  f0

coefficient = -j*2*pi;

for i= 1:P
    tao(i) = (i-1)*d*sind(theta)/c;
end

a = exp(coefficient*f0*tao);

Un = V(:,M+1:P);

J = a*Un*(a*Un)';
J = 1/abs(J);