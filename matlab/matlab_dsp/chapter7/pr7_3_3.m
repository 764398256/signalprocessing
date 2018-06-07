% 
% pr7_3_3  
clear all; clc; close all;

f0=49.13;                    % ����Ƶ��
fs=3200;                     % ����Ƶ��
N=2048;                      % ���ݳ���
n=0:N-1;                     % ��������
rad=pi/180;                  % �ǶȺͻ��ȵ�ת������
xb=[240,0.1,12,0.1,2.7,0.05,2.1,0,0.3,0,0.6]; % г����ֵ
Q=[0,10,20,30,40,50,60,0,80,0,100]*rad;       % г����ʼ��λ

M=11;                        % г������
df=fs/N;                     % Ƶ�ʷְ���
t=n/fs;                      % ʱ������
% Blackman-Harris��
win=blackmanharris(N)';
x=zeros(1,N);                % ��ʼ��
for k=1 : M                  % ����г���ź�
x=x+xb(k)*cos(2*pi*f0*k*t+Q(k));
end
X=fft(x.*win);               % �źų��Դ�������FFT
Y=abs(X);                    % ȡƵ�׵ķ�ֵ
A=zeros(1,M);                % ��ʼ��
ff=zeros(1,M);
Ph=zeros(1,M);

for k=1 : M                  % ��������͸���г���Ĳ���
    if k==1                  % ���������,��40-60Hz������Ѱ������ֵ
        n1=fix(40/df); n2=fix(60/df); % ���40Hz��60Hz��Ӧ��������
    else                     % ������г��,�Ӹ�г������ֵ-10��+10��������Ѱ�����ֵ
        n1=fix((k*ff(1)-10)/df); % ��������Ӧ��������
        n2=fix((k*ff(1)+10)/df);
    end
    [Fm,nn]=max(Y(n1:n2));   % ���������ҳ����ֵ
    nm=nn+n1-1;              % �������ֵ��������
    if Y(nm+1)==Y(nm-1)
        delta=0;
    elseif Y(nm+1)<Y(nm-1);
        nm=nm-1;
    end
    y1=Y(nm);                % ���y1��y2
    y2=Y(nm+1);
    beta=(y2-y1)/(y2+y1);    % ��ʽ(7-3-3)�����beta
% Blackman-Harris����beta��alpha�ı�ʾʽ
    alpha=2.61979588*beta+0.286567*beta^3+0.12830543*beta^5+0.0802152*beta^7;
% ��ʽ(7-3-9)��Blackman-Harris���Ĵ�������ϵ�����г���ķ�ֵ
    A(k)=(y1+y2)*(3.06539914+0.96556547*alpha^2+0.163418995*alpha^4+...
       0.02080189*alpha^6)/2;
    ff(k)=(nm-1+alpha+0.5)*fs/N;         % ���г����Ƶ��
    Ph(k)=angle(X(nm))-pi*(alpha+0.5);   % ���г���ĳ�ʼ��λ
    Ph(k)=Ph(k)-(Ph(k)>pi)*2*pi+(Ph(k)<-pi)*2*pi; % ����λ��������
% ����ֵ��С��Ϊ0,����Ƶ����λ����
    if A(k)<0.0005, A(k)=0; ff(k)=k*ff(1); Ph(k)=0; end    
% ��ʾг������
    fprintf('%4d      %5.6f  %5.6f  %5.6f\n',k,ff(k),A(k),Ph(k)/rad);
end
