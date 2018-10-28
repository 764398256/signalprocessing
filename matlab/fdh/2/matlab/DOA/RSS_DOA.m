function J = RSS_DOA(V,f0,theta)
% Circle array music J function
% theta
% N              the length of data
% R              circle array senser number
% Fs             signal frequency
% r              circle array's radius
% c              velocity of sound
% L              the number of samples in a segment
global  P   d  c M 

coefficient = j*2*pi;

J = 0;

for i = 1:P
%     tao(i) = r*cos(2*pi*(i-1)/R-theta)/c;
    tao(i) = (i-1)*d*sin(theta)/c;
    if i<=P-M
        G(:,i) = V(:,i+M);
    end
end
a = exp(coefficient*f0*tao);
J = J+abs(a*G*G'*a');
J = 1/J;



