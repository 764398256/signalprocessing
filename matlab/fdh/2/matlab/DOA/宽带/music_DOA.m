function J = music_DOA(V,theta)

global  P   d  c M Fl Fh

coefficient = j*2*pi;

J = 0;
bin_num = size(V,3);
deltf = (Fh-Fl)/(bin_num-1);
for n = 1:bin_num
    for i = 1:P
        tao(i) = (i-1)*d*sin(theta)/c;
        if i<=P-M
            G(:,i) = V(:,i+M,n);
        end
    end
    f = Fl+deltf*(n-1);
    a = exp(coefficient*f*tao);
    J = J+abs(a*G*G'*a');
end
J = 1/J;



