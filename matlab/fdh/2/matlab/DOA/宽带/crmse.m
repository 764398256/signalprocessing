function rmse = crmse(peak,theta)

M = length(theta);
theta = theta*180/pi;
rmse = 0;
rmse_temp = zeros(1,M);

for m = 1:M
    for i = 1:M
       rmse_temp(i) = abs(theta(m)-peak(1,i))^2;
    end
    rmse = rmse+min(min(rmse_temp));
end