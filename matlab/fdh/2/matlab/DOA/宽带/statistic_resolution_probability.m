clear all

global  P   d  c M Fl Fh

c = 340;
P = 8;
M = 2;
bin_num = 9;
snapshot = 50;
Fl = 80;
Fh = 120;
f0 = 100;
d = c/(2*f0);
amplitude = [1,1];
theta0 = [11,13]*pi/180;

tic
Pro = zeros(3,20);
p_num = 0;
for SNR = -10:5:30
    SNR
    km = 0;
    kr = 0;
    kt = 0;
    EXP_NUM = 200;
    for exp_num = 1:EXP_NUM
        X = generate_signal_Wbins(Fl,Fh,bin_num,snapshot,SNR,amplitude,theta0);

        theta = [9,12,15]*pi/180;

%         [V,D] = music_eig(X);

        [V_r,D_r] = RSS_eig(X,theta,f0);

        [V_T,D_T,f0_T] = TCT_eig(X,theta);

        step = 0.01;
        Jm = zeros(1,180/step);
        Jr = Jm;
        Jt = Jm;

        t = 0;
        for k = 1:step:180.99
            t = t+1;
            theta = (k-90)*pi/180;

%             Jm(t) = music_DOA(V,theta);
            Jr(t) = RSS_DOA(V_r,f0,theta);
            Jt(t) = TCT_DOA(V_T,f0_T,theta);

        end
        delt = abs(theta0(1)-theta0(2))/2;
        delt = delt*180/pi;

        peak_Jm = peak_find_rmse(Jm,step);
        delt1_m = abs(peak_Jm(1,1)-theta0(1)*180/pi);
        delt2_m = abs(peak_Jm(1,2)-theta0(2)*180/pi);
        if delt1_m < delt && delt2_m < delt && peak_Jm(1,2)~=0
            km = km+1;
        end

        peak_Jr = peak_find_rmse(Jr,step);
        delt1_r = abs(peak_Jr(1,1)-theta0(1)*180/pi);
        delt2_r = abs(peak_Jr(1,2)-theta0(2)*180/pi);
        if delt1_r < delt && delt2_r < delt && peak_Jr(1,2)~=0
            kr = kr+1;
        end

        peak_Jt = peak_find_rmse(Jt,step);
        delt1_t = abs(peak_Jt(1,1)-theta0(1)*180/pi);
        delt2_t = abs(peak_Jt(1,2)-theta0(2)*180/pi);
        if delt1_t < delt && delt2_t < delt && peak_Jt(1,2)~=0
            kt = kt+1;
        end
    end
    p_num = p_num+1;
    Pro(1,p_num) = km/EXP_NUM;
    Pro(2,p_num) = kr/EXP_NUM;
    Pro(3,p_num) = kt/EXP_NUM;
end
toc
