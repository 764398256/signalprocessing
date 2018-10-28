clear all

data = wavread('sai_4mic.wav');
fs = 16000;
MicNum = 4;
MArray = [0.035,0;0,0.035;-0.035,0;0,-0.035];

N_L = 4096;
N = 256;
Snap = N_L/N;
n_blocks = floor(length(data)/N_L);

angle_beam = zeros(n_blocks,1);
angle_music = zeros(n_blocks,1);
for i = 1:n_blocks
    N1 = (i-1)*N_L+1;
    N2 = i*N_L;
    x = data(N1:N2,:);
    
    Rx = zeros(MicNum,MicNum,N/2);
    for n = 1:Snap
        X = fft(x((n-1)*N+1:n*N,:));
        for k = 1:N/2
            Rx(:,:,k) = Rx(:,:,k)+X(k,:)'*X(k,:);
        end
    end
    Rx = Rx./Snap;
    
    V = zeros(MicNum,MicNum,N/2);
    for k = 1:N/2
        [V(:,:,k),D] = eig_sort(Rx(:,:,k));
    end
    
    J_beam = zeros(360,1);
    J_music = zeros(360,1);
    for theta = 1:360
        J_beam(theta) = wideband_beamform_doa_real(Rx,theta,MArray);
        J_music(theta) = wideband_music_doa_real(V,theta,MArray);
    end
    
    [val,I] = max(J_beam);
    angle_beam(i) = I;
    [val,I] = max(J_music);
    angle_music(i) = I;
end
[B, f, t] = spectrogram(data(:,1)*10,1024,512,1024,fs);
imagesc(t,f,20*log10(abs(B))); 
axis xy;

figure
t = (1:n_blocks)*4096/fs;
plot(t,angle_beam);
hold on
plot(t,angle_music,'r');
    
    
        