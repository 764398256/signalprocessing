function peak = peak_find(input,step,search_range)

Threshold = 1;
search_begin = search_range+1;
search_end = length(input)-search_range;

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
% k = 1;
% for i = search_begin:search_end
%     if (input(i)>input(i-1) && input(i)>input(i+1))
%         delt1 = input(i)-min(input(i-search_range : i-1));
%         delt2 = input(i)-min(input(i+1 : i+search_range));
%         if delt1 > Threshold && delt2 > Threshold 
%             theta = i*step-(89+step);
%             peak(1,k) = theta;
%             peak(2,k) = input(i);
%             k = k+1;
%         end
%     end
% end

