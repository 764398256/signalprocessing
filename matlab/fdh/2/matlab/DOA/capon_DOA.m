function J = capon_DOA(Rx_1,theta)

global  P   d  c M  f0

coefficient = -j*2*pi;

for i= 1:P
    tao(i) = (i-1)*d*sind(theta)/c;
end

a = exp(coefficient*f0*tao);

J = abs(a*Rx_1*a');
J = 1/J;