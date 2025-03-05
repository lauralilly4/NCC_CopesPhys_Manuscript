%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOAA/OSU post-doc: Figs. 6, S1 %%%
%%%     - ALF transitions          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laura Lilly (code adapted from B.T. Cervantes)
% Updated: 23 Jul 2023
% Plot: - yearly plots of daily-averaged ALF at NH10 -> to determine
% transition date of 'Raw' ALF 


%% ======= Load velocity file & assign variables =======  
% loadfile = 'Phys_Inds/BTC_ADCP_NH10_1997_2021_V6.mat';
loadfile = 'Phys_files_2024/ADCP_NH10_Rotated_1997_2024_V5.mat';
load(loadfile)
clear sal temp

% Depth-average velocities
udavg = nanmean(u,1);
vdavg = nanmean(v,1);
clear u v
% Low-pass filter depth-averaged velocities
window = 36;
udavg_filt = pl66tn(fillmissing(udavg,'linear'),1,window); 
vdavg_filt = pl66tn(fillmissing(vdavg,'linear'),1,window);
udavg_filt(find(isnan(udavg))) = NaN;
vdavg_filt(find(isnan(vdavg))) = NaN;
clear window
% Use low-pass filtered velocities to calculate principal axes
a = udavg_filt;  
b = vdavg_filt;  
a(isnan(b))=nan;
b(isnan(a))=nan;
a=a(~isnan(a));
b=b(~isnan(b));
z = complex(a,b);
[theta,~,~,~]=princax(z);
[along,across] = rotate2princax(udavg,vdavg,theta*180/pi);
[along_filt,across_filt] = rotate2princax(udavg_filt,vdavg_filt,theta*180/pi);
clear a b z theta 
% Calculate daily mean, standard deviation, and confidence intervals
conf = 0.95;
[along_N,along_mean,along_std,along_ci] = daily_clim(conf,along,time);
% Calculate 2-week filtered velocities
window = 24*14;
along_2wfilt = pl66tn(fillmissing(along,'linear'),1,window); 
across_2wfilt = pl66tn(fillmissing(across,'linear'),1,window);
along_2wfilt(find(isnan(along))) = NaN;
across_2wfilt(find(isnan(across))) = NaN;
clear window

% % Plot Annual Along-shelf Velocities with climatology and 95% confidence
% % interval
% figure('position',[25 450 1000 1500])
% ref_yrday = datenum(1998,1,1):1:datenum(1998,12,31);


%% ======= PLOT: Single year-chunk of ALF vlaues =======
yrin = input('Which year?   ');
yeardnum = datenum(yrin,1,1);
% pst_ind = find(pst_dts.Year == yrin);
% bst_ind = find(bst_dts.Year == yrin);
ref_yrday = datenum(yrin,1,1):1:datenum(yrin,12,31);
ref_dts = datetime(ref_yrday,'ConvertFrom','datenum');

if sum(yrin == [2000,2004,2008,2012,2016,2020]) > 0
    ref_dts(datefind(datetime(yrin,2,29,0,0,0),ref_dts)) = [];
end


figure(1)
p = ciplot(along_mean-along_std,along_mean+along_std,ref_yrday(1:365),[225/255 225/255 252/255]); 
hold on; 
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
plot(ref_yrday(1:365),along_mean,'color',[200/255 200/255 252/255],'linewidth',1); 
ref_year = (str2num(datestr(yeardnum,10))-str2num(datestr(ref_yrday(1:365),10)))*365.25; 
ind = find(time >= datenum(yrin,1,1) & time < datenum(yrin+1,1,1));
plot(time(ind)-ref_year,along_filt(ind),'b','linewidth',1);
plot(time(ind)-ref_year,along_2wfilt(ind),'b','linewidth',2);
plot(time(ind)-ref_year,zeros(1,length(ind)),'k');

set(gca,'ylim',[-0.8 0.8]);
x_start = time(ind(1))-ref_year(1); 
x_end = time(ind(end))-ref_year(end);
set(gca,'xlim',[x_start x_end]);
set(gca,'xtick',[x_start x_start+32 x_start+60 x_start+91 x_start+121 ...
    x_start+152 x_start+182 x_start+213 x_start+244 x_start+274 x_start+305 ...
    x_start+335])
set(gca,'xticklabel',['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep'; ...
    'Oct';'Nov';'Dec'])
set(gca,'fontsize',12);
ylabel('Alongshore Flow Mag (m/s)','fontsize',12);
set(gca,'box','off','TickDir','out');
set(gcf,'color','white');


keyboard


% %%% Addendum #1: Get *daily avgs* of instantaneous flow -> to align with
% cumulative flows
flwtim_dt = datetime(time','ConvertFrom','datenum'); 
instflw_hr = timetable(flwtim_dt, along_filt);
instflw_dy = retime(instflw_hr, 'daily', @(x) mean(x, 'omitnan'));



% %%% Addendum #2: Table of 'along_mag' and 'date' -> to ID beginning of
% sustained southward flow
dts_mag_tbl = table(ref_dts',along_mean);

filt_tim = time(ind)-ref_year;
filt_dts = datetime(filt_tim(1,:),'ConvertFrom','datenum');
magfilt = along_2wfilt(ind);
dtname = input('Which Mon-Dy?   [mm,dd]  ');
dtst = datetime(yrin,dtname(1),dtname(2),0,0,0);

dtid = datefind(dtst,filt_dts);

mag_filt_tbl = table(filt_dts(dtid:end)',magfilt(dtid:end));

mag_filt_tbl(1:1000,:)
