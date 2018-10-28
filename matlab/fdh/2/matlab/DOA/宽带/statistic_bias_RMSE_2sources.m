clear all

global  P   d  c M Fl Fh

c = 340;
P = 8;
M = 2;
bin_num = 33;
snapshot = 100;
Fl = 80;
Fh = 120;
f0 = 100;
d = c/(2*f0);
amplitude = [1,1];

tic
bias_ave = zeros(3,50,10);
rmse_ave = zeros(3,50,10);

theta0 = [11,13]*pi/180;
sigma_theta = [1];
for sigma_sequence = 1:length(sigma_theta)
    p_num = 0;
    sigma_sequence
    for SNR = 5:5:25
        SNR
        rmse_m = 0;

        rmse_r = 0;

        rmse_t = 0;
        
        EXP_NUM = 200;
        for exp_num = 1:EXP_NUM
            X = generate_signal_Wbins(Fl,Fh,bin_num,snapshot,SNR,amplitude,theta0);

            theta = [9,12,15]*pi/180;

            [V,D] = music_eig(X);

            [V_r,D_r] = RSS_eig(X,theta,f0);

%             [V_T,D_T,f0_T] = TCT_eig(X,theta);

            step = 0.01;
            Jm = zeros(1,180/step);
            Jr = Jm;
            Jt = Jm;

            t = 0;
            for k = 1:step:180.99
                t = t+1;
                theta = (k-90)*pi/180;

                Jm(t) = music_DOA(V,theta);
                Jr(t) = RSS_DOA(V_r,f0,theta);
%                 Jt(t) = TCT_DOA(V_T,f0_T,theta);

            end
            
            peak_Jm = peak_find_rmse(Jm,step);
            rmse_m_temp = crmse(peak_Jm,theta0);
            rmse_m = rmse_m+rmse_m_temp;

            peak_Jr = peak_find_rmse(Jr,step);
            rmse_r_temp = crmse(peak_Jr,theta0);
            rmse_r = rmse_r+rmse_r_temp;
            
%             peak_Jt = peak_find_rmse(Jt,step);
%             rmse_t_temp = crmse(peak_Jt,theta0);
%             rmse_t = rmse_t+rmse_t_temp;


        end
        p_num = p_num+1;

        rmse_ave(1,p_num,sigma_sequence) = sqrt(rmse_m/EXP_NUM/M);
        rmse_ave(2,p_num,sigma_sequence) = sqrt(rmse_r/EXP_NUM/M);
        rmse_ave(3,p_num,sigma_sequence) = sqrt(rmse_t/EXP_NUM/M);
    end
end
toc
