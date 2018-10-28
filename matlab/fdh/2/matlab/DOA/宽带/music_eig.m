function [V,D] = music_eig(X)
% Circle array music J function
% theta
% N              the length of data
% R              circle array senser number
% Fs             signal frequency
% r              circle array's radius
% c              velocity of sound
% L              the number of samples in a segment
global  P

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