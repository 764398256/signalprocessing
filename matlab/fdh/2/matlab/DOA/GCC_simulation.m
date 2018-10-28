clear all

data = wavread('4mic_simulate.wav');
fs = 16000;
MicNum = 4;
MArray = [0.1,0;0,0.1;-0.1,0;0,-0.1];

% data = data+0.001*randn(length(data),MicNum);

N_L = 1024;
N = 256;
Snap = N_L/N;
n_blocks = floor(length(data)/N_L);
c = 340;

angle_gcc = zeros(n_blocks,1);
for i = 1:n_blocks
    N1 = (i-1)*N_L+1;
    N2 = i*N_L;
    x = data(N1:N2,:);
    X = fft(x);
    if i == 15
        i
    end
    for chm = 1:MicNum-1
        for chn = chm+1:MicNum
            Rmn = X(:,chm).*conj(X(:,chn));
%             Rt = ifft(Rmn); %GCC
            Rt = ifft(Rmn./abs(Rmn)); %PHAT
            plot(Rt);
            
            step = 0.01;
            xx=[1:step:N_L];
            x=[1:N_L];
            Rtp=spline(x,Rt,xx);
            
            [peak,peak_id] = max(Rtp);
            peak_id = (peak_id - step)*step;
            if peak_id < N_L/2
                delay(chm,chn) = peak_id;
            end
            if peak_id > N_L/2
                delay(chm,chn) = -(N_L-peak_id);
            end
        end
    end
    
    delay = delay/fs*c;
    dMArray = zeros(MicNum*(MicNum-1)/2,2);
    dDelay = zeros(MicNum*(MicNum-1)/2,1);
    
    n = 1;
    for chm = 1:MicNum-1
        for chn = chm+1:MicNum
            dMArray(n,:) = (MArray(chm,:) - MArray(chn,:));
            dDelay(n) = -delay(chm,chn);
            n = n + 1;
        end
    end
    
    step = 1;
    J1 = zeros(360/step,1);
    Jmin = 1000;
    k = 0;
    for phi = 1:step:360
        k = k+1;
        phi1 = phi*pi/180;
        J1(k) = sum((abs(dDelay-dMArray*[cos(phi1),sin(phi1)].')).^2);
        J = J1(k);
        if J < Jmin
            Jmin = J;
            phi_f = phi1;
        end
    end

    angle_gcc(i) = phi_f*180/pi;
end
[B, f, t] = spectrogram(data(:,1)*10,1024,512,1024,fs);
imagesc(t,f,20*log10(abs(B)));
axis xy;

figure
t = (1:n_blocks)*N_L/fs;
plot(t,angle_gcc,'o');


