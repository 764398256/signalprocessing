function J = conventional_DOA(Rx,theta)

global  P   d  c M  f0

coefficient = -j*2*pi;

for i= 1:P
    tao(i) = (i-1)*d*sind(theta)/c;
end

a = exp(coefficient*f0*tao);

J = abs(a*Rx*a');