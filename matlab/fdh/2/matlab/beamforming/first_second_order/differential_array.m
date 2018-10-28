clc;
clear all;
close all;
%%
x_all = audioread('sai_3mic.wav');
fs =16000;
FrmLen = 512;
FrmShift = 256;
AnaWnd = hanning(FrmLen);
AnaWnd = sqrt(AnaWnd);
FreBin = (0 : fs/FrmLen : fs/2)';

d = 0.035;% m
c = 340; % m/s
alpha = 1/5; % alpha =0,dipole;1/2 cardioid
order = 2;% 1:first order; 2:second order
%%
% for special case : uniformly-spaced microphones array
%               |
%               |
%           d   |   d
%--------1------2------3------->
%               |
%               |
%               |
%% plot
theta = 0:2*pi/360:2*pi-2*pi/360;
bf = alpha + (1-alpha)*cos(theta); %first order
polar(theta,abs(bf),'r');
hold on
bf = (alpha + (1-alpha)*cos(theta)).^2; %second order
polar(theta,abs(bf));
%% first order,and adopt 1&3 mic
if order == 1
    d_1_2 = d;
    tau = d_1_2/c;
    tau_type = alpha*tau; 
    kd = 2 * pi * FreBin * (tau_type);
    comp = abs(1./(1-exp(-j*2*pi*FreBin*(tau_type+tau))));
    comp(1) = comp(2);
    comp = [comp;flipud(comp(2:end-1,:))];
    comp = min(comp,20);
%     comp = ones(FrmLen,1);
    x = [x_all(:,1),x_all(:,2)];
    clear x_all
    len = length(x);
    xinp = zeros(FrmLen,2);
    direct1 = zeros(len,1);
%     direct3= zeros(len,1);
    FrmFlag = 1;
    TimeIdx = 1;
    while TimeIdx + FrmShift - 1 <= len
        xinp = [xinp(FrmShift+1:end,:);x(TimeIdx:TimeIdx+FrmShift-1,1:2)];

        Xfft = fft(xinp.*repmat(AnaWnd,1,2));

        %% Delay version of the input signal
        tmpYfft = Xfft(1:FrmLen/2+1,:).*repmat(exp(-1j*kd),1,2);
        Yfft  = [tmpYfft;flipud(conj(tmpYfft(2:end-1,:)))];
        
        tmp = Xfft(:,1) - Yfft(:,2);
        tmp = tmp.*comp;
        zfrm(:,1) = real(ifft(tmp));

        if(FrmFlag == 1),
            direct1(TimeIdx : 1 : TimeIdx + FrmShift - 1,1) = zfrm(1 : FrmShift,1);
        else
            direct1(TimeIdx : 1 : TimeIdx + FrmShift - 1,1) = Prezfrm(end - FrmShift + 1 : end,1) +  zfrm(1 : FrmShift,1);
        end;
        Prezfrm = zfrm;
        FrmFlag = FrmFlag + 1;
        TimeIdx = TimeIdx + FrmShift;
    end;
     direct1 = direct1(FrmShift+1:end);
    wavwrite(direct1./max(abs(direct1)),fs,16,'Carioid_first_c.wav');
%%
% second order,and adopt 1&5&3 mic
else
    d_2_3 = d;
    d_1_2 = d;
    tau = d/c;
    tau_type = alpha*tau; 
    kd = 2 * pi * FreBin * (tau_type);
    comp = abs(1./(1-exp(-j*2*pi*FreBin*(tau_type+tau))));
    comp(1) = comp(2);
    comp = comp.^2;
    comp = [comp;flipud(comp(2:end-1,:))];
    comp = min(comp,20);
%     comp = ones(FrmLen,1);
    x = [x_all(:,1),x_all(:,2),x_all(:,3)];
    clear x_all
    len = length(x);
    xinp = zeros(FrmLen,3);
    direct1 = zeros(len,1);
    direct3 = zeros(len,1);
    FrmFlag = 1;
    TimeIdx = 1;
    while TimeIdx + FrmShift - 1 <= len
        xinp = [xinp(FrmShift+1:end,:);x(TimeIdx:TimeIdx+FrmShift-1,1:3)];

        Xfft = fft(xinp.*repmat(AnaWnd,1,3));

        %% Delay version of the input signal
        tmpYfft = Xfft(1:FrmLen/2+1,1:3).*repmat(exp(-1j*kd),1,3);
        Yfft  = [tmpYfft;flipud(conj(tmpYfft(2:end-1,:)))];
        
        tmp = (Xfft(:,1) -Yfft(:,2) - (Xfft(:,2)-Yfft(:,3)));
        tmp = tmp.*comp;
        zfrm(:,1) = real(ifft(tmp));

        if(FrmFlag == 1),
            direct1(TimeIdx : 1 : TimeIdx + FrmShift - 1,1) = zfrm(1 : FrmShift,1);
        else
            direct1(TimeIdx : 1 : TimeIdx + FrmShift - 1,1) = Prezfrm(end - FrmShift + 1 : end,1) +  zfrm(1 : FrmShift,1);
        end;
        Prezfrm = zfrm;
        FrmFlag = FrmFlag + 1;
        TimeIdx = TimeIdx + FrmShift;
    end;
    direct1 = direct1(FrmShift+1:end);
    wavwrite(direct1./max(abs(direct1)),fs,16,'Carioid_second_c.wav');
end
    





