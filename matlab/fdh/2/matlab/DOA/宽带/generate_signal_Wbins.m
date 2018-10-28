function X = generate_signal_Wbins(Fl,Fh,bin_num,snapshot,SNR,amplitude,theta)
% N              FFT length
% R              circle array senser number
% Fs             signal frequency
% r              circle array's radius
% c              velocity of sound
global  P   d  c M 

f = Fl:(Fh-Fl)/(bin_num-1):Fh;
coefficient = -j*2*pi;
X = zeros(P,snapshot,length(f));
for k = 1:length(f)
    S = zeros(M,snapshot);
    for i = 1:M
        S(i,:) = amplitude(i)*(randn(1,snapshot)+j*randn(1,snapshot));
    end
    for i = 1:P
        for m = 1:M
            tao(i,m) = (i-1)*d*sin(theta(m))/c;
        end
    end
    sigma = amplitude(1)/(10^(SNR/20));
    A = exp(coefficient*f(k)*tao);
    for n = 1:snapshot
        W = sigma*(randn(P,1)+j*randn(P,1));
        X(:,n,k) = A*S(:,n)+W;
    end
end

