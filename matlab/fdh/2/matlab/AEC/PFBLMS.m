clear all

%
load('wo.mat');
w0 = wo;
plot(w0)
h = [0.35,1,-0.35];

runs = 20;
sigman = 0.1;
itn = 300000/2;
N = length(w0);

xi=zeros(itn,1);
M = 64;
P = N/M;
mu = 0.5;
alpha = 0.99;
for k = 1:runs
	x = filter(h,1,randn(itn,1));
    tmp = conv(x,w0);
	d = tmp(1:itn)+sigman*randn(itn,1);
    wF = zeros(2*M,P);
    xF = zeros(2*M,P);
    block_n = 1;
    for n = M+1:M:length(x)-M
        xF = [fft(x(n-M:n+M-1)),xF(:,[1:end-1])];
        yhat = ifft(sum((wF.*xF).').');
        yhat = real(yhat(M+1:end));
        E = d(n:n+M-1)-yhat;
        MU = mu./(sum((abs(xF).^2)')'+0.1);
        EF = fft([zeros(M,1);E]);
        for p = 1:P
            wF(:,p) = wF(:,p)+MU.*EF.*conj(xF(:,p));
        end
        waux = real(ifft(wF));
        wF = fft([waux(1:M,:); zeros(M,P)]);
        xi(block_n) = xi(block_n)+sum(E.^2)/M;
        block_n = block_n+1;
    end
end
xi=xi/runs;
semilogy(xi)

xi=zeros(itn,1);
M = 128;
P = N/M;
mu = 0.5;
alpha = 0.99;
for k = 1:runs
	x = filter(h,1,randn(itn,1));
    tmp = conv(x,w0);
	d = tmp(1:itn)+sigman*randn(itn,1);
    wF = zeros(2*M,P);
    xF = zeros(2*M,P);
    block_n = 1;
    for n = M+1:M:length(x)-M
        xF = [fft(x(n-M:n+M-1)),xF(:,[1:end-1])];
        yhat = ifft(sum((wF.*xF).').');
        yhat = real(yhat(M+1:end));
        E = d(n:n+M-1)-yhat;
        MU = mu./(sum((abs(xF).^2)')'+0.1);
        EF = fft([zeros(M,1);E]);
        for p = 1:P
            wF(:,p) = wF(:,p)+MU.*EF.*conj(xF(:,p));
        end
        waux = real(ifft(wF));
        wF = fft([waux(1:M,:); zeros(M,P)]);
        xi(block_n) = xi(block_n)+sum(E.^2)/M;
        block_n = block_n+1;
    end
end
xi=xi/runs;
hold on
semilogy(xi,'g')

xi=zeros(itn,1);
M = 256;
P = N/M;
mu = 0.5;
alpha = 0.99;
for k = 1:runs
	x = filter(h,1,randn(itn,1));
    tmp = conv(x,w0);
	d = tmp(1:itn)+sigman*randn(itn,1);
    wF = zeros(2*M,P);
    xF = zeros(2*M,P);
    block_n = 1;
    for n = M+1:M:length(x)-M
        xF = [fft(x(n-M:n+M-1)),xF(:,[1:end-1])];
        yhat = ifft(sum((wF.*xF).').');
        yhat = real(yhat(M+1:end));
        E = d(n:n+M-1)-yhat;
        MU = mu./(sum((abs(xF).^2)')'+0.1);
        EF = fft([zeros(M,1);E]);
        for p = 1:P
            wF(:,p) = wF(:,p)+MU.*EF.*conj(xF(:,p));
        end
        waux = real(ifft(wF));
        wF = fft([waux(1:M,:); zeros(M,P)]);
        xi(block_n) = xi(block_n)+sum(E.^2)/M;
        block_n = block_n+1;
    end
end
xi=xi/runs;
hold on
semilogy(xi,'r')

