clear alll

c = 340;
d = 0.17;
P = 9;
MArray = [-4*d,-3*d,-2*d,-d,0,d,2*d,3*d,4*d]';

theta = 0;
A = exp(-j*2*pi*f*MArray*sind(theta)/c);
win = ones(length(A),1)/length(A);
Aw = A.*win;
J = zeros(181,50);

f_idx = 1;
for f = 50:50:2500
    idx = 1;
    for theta = -90:1:90
        a = exp(-j*2*pi*f*MArray*sind(theta)/c);
        J(idx,f_idx) = abs(Aw'*a);
        idx = idx+1;
    end
    plot((J(:,f_idx)))
    title(['frequency:  ',num2str(f)]);
    pause(0.2)
    clc
    f_idx = f_idx+1;
end