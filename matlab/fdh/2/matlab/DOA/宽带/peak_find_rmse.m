function peak = peak_find_rmse(input,step)


input = 10*log10(input./max(max(input)));
peak = 200*ones(2,10);

[peak_val,peak_loc] = findpeaks(input,'sortstr','descend');
for i = 1:length(peak_val)
    peak(1,i) = peak_loc(i)*step-(89+step);
    peak(2,i) = peak_val(i);
end

if peak(1,1) > peak(1,2)
    peak_temp = peak(:,1);
    peak(:,1) = peak(:,2);
    peak(:,2) = peak_temp;
end


