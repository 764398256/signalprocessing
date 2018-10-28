clear all

%
w0 = ones(15,1);
w0(9:15) = -1*w0(9:15);

h = [0.35,1,-0.35];

runs = 100;
sigman = 0.1;
itn = 3000;
xi=zeros(itn,1);
N = length(w0);

mu = 0.5;
alpha = 0.99;
for k=1:runs
	x=filter(h,1,randn(itn,1));
	d=filter(w0,1,x)+sigman*randn(itn,1);
	w=zeros(N,1);
    xtdl=zeros(size(w));
    pxx = 1;
	for n=1:itn
		xtdl=[x(n);xtdl(1:length(xtdl)-1)];
		yhat=w'*xtdl;
		e=d(n)-yhat;
%         pxx = xtdl'*xtdl;
%         w=w+mu*e*xtdl/(pxx+0.01);
        pxx = alpha*pxx+(1-alpha)*xtdl(1)*xtdl(1);
		w=w+mu*e*xtdl/(N*pxx+0.01);
		xi(n)=xi(n)+e^2;
	end
end
xi=xi/runs;
semilogy(xi)

mu = 0.02;
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
hold on
semilogy(xi,'r')


delta = 0.0001;
lambda = 0.985;

xi=zeros(itn,1);
for k=1:runs
	x=filter(h,1,randn(itn,1));
	d=filter(w0,1,x)+sigman*randn(itn,1);
	w=zeros(N,1);
    xtdl=zeros(size(w));
	Rxx_inv=(1/delta)*eye(N);
	for n=1:itn
		xtdl=[x(n);xtdl(1:length(xtdl)-1)];
		u=Rxx_inv*xtdl;
		k=u/(lambda+xtdl'*u);
		yhat=w'*xtdl;
		e=d(n)-yhat;
		w=w+k*e;
		Rxx_inv=(1/lambda)*(Rxx_inv-k*(xtdl'*Rxx_inv));
		xi(n)=xi(n)+e^2;
	end
end
xi=xi/runs;
semilogy(xi,'g')

