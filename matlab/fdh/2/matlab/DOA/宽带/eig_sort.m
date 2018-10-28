function [V,D] = eig_sort(RX)
% V:the eigenvectors of RX
% D:the eigenvalues of RX in descending order 

[V1,D1] = eig(RX);
D2 = ones(1,length(D1))*D1;
[D,I] = sort(D2,'descend');

% [V1,D1] = jacob_H(RX);
% [D,I] = sort(D1,'descend');
for i = 1:length(D1)
    V(:,i) = V1(:,I(i));
end