function [N,daily_mean,stddev,ci] = daily_clim(conf,data,time)

% Calculate means for each year day for an annual climatology 
% and confidence intervals on mean values
ref_yrday = datenum(2000,01,01):1:datenum(2000,12,31); 
time_yrday = cellstr(datestr(time,6)); 
N = NaN(length(ref_yrday),1); 
daily_mean = NaN(length(ref_yrday),1); 
stddev = NaN(length(ref_yrday),1); 
ci = NaN(length(ref_yrday),2); 
alpha = 1-conf; pLo = alpha/2; pUp = 1-alpha/2;
    for i = 1:366
        ind = find(strcmp(time_yrday,datestr(ref_yrday(i),6)));
        good = find(~isnan(data(ind)));
        N(i) = length(unique(str2num(datestr(time(ind(good)),10))));
        if ~isempty(ind)
            daily_mean(i) = nanmean(data(ind));
            stddev(i) = nanstd(data(ind));
            se = nanstd(data(ind))/sqrt(N(i));
            crit = tinv([pLo pUp],N(i)-1);
            ci(i,:) = daily_mean(i) + crit*se;
        end
    end
N = [N(1:59); N(61:366)];
daily_mean = [daily_mean(1:59); daily_mean(61:366)];
stddev = [stddev(1:59); stddev(61:366)];
ci = [ci(1:59,:); ci(61:366,:)];