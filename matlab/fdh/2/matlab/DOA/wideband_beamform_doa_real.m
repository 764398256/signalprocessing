function J = wideband_beamform_doa_real(Rx,theta,MArray)

fs = 16000;
N = 256;
c = 340;

coefficient = -j*2*pi*fs/(N*c);
h = [cosd(theta),sind(theta)].';
tao = MArray*h;

D = exp(coefficient*tao);
J = 0;
for k = 1:N/2
    A = D.^(k-1);
    J = J+abs(A'*Rx(:,:,k)*A);
end