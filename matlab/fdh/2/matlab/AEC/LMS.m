clear all

%´ý¹À¼ÆÂË²¨Æ÷
w0 = ones(15,1);
w0(9:15) = -1*w0(9:15);
% freqz(w0,1)

%ÓÐÉ«ÔëÉùÂË²¨Æ÷
h = [0.35,1,-0.35];
% freqz(h,1)

runs = 100;
sigman = 0.1;
itn = 3000;
xi = zeros(itn,1);
N = length(w0);
mu = 0.002;
for k=1:runs
    x=filter(h,1,randn(itn,1));
    d=filter(w0,1,x)+sigman*randn(itn,1);
    w=zeros(N,1);
    xtdl=zeros(size(w));
    for n=1:itn
        xtdl=[x(n);xtdl(1:length(xtdl)-1)];
        yhat=w'*xtdl;
        e=d(n)-yhat;
        w=w+2*mu*e*xtdl;
        xi(n)=xi(n)+e^2;
    end
end
xi=xi/runs;
semilogy(xi)
hold on

mu = 0.01;
for k=1:runs
    x=filter(h,1,randn(itn,1));
    d=filter(w0,1,x)+sigman*randn(itn,1);
    w=zeros(N,1);
    xtdl=zeros(size(w));
    for n=1:itn
        xtdl=[x(n);xtdl(1:length(xtdl)-1)];
        yhat=w'*xtdl;
        e=d(n)-yhat;
        w=w+2*mu*e*xtdl;
        xi(n)=xi(n)+e^2;
    end
end
xi=xi/runs;
semilogy(xi,'r')

mu = 0.03;
for k=1:runs
    x=filter(h,1,randn(itn,1));
    d=filter(w0,1,x)+sigman*randn(itn,1);
    w=zeros(N,1);
    xtdl=zeros(size(w));
    for n=1:itn
        xtdl=[x(n);xtdl(1:length(xtdl)-1)];
        yhat=w'*xtdl;
        e=d(n)-yhat;
        w=w+2*mu*e*xtdl;
        xi(n)=xi(n)+e^2;
    end
end
xi=xi/runs;
semilogy(xi,'g')

