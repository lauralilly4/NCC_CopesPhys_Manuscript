%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOAA/OSU post-doc: Fig. 8b,d (new)  %%%
%%%       - ALF transitions             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laura Lilly (code adapted from B.T. Cervantes)
% Updated: 14 Jul 2024
% Extract and combine 'instantaneous' and 'cumulative' alongshore flows at
% NH10, using a common daily timeseries (1997-2023)


%% ======== File load ========
% %%% Instantaneous flow file
instfl = 'Phys_files_2024/ADCP_NH10_Rotated_1997_2024_V5.mat';
load(instfl)
along_inst = nanmean(along,1);
acrs_inst = nanmean(across,1);
time_inst = time;
dt_inst = datetime(time,'ConvertFrom','datenum');

% %%% Cumulative flow file
cumufl = load('Phys_files_2024/Wind_Indices_Vel_45N_V2_1JanCumFlow.mat');
flow_cumu = cumufl.cum_vel;
time_cumu = cumufl.time_vel;


%% ======== Adjust and combine variables ========
% First, convert instantaneous flow interval from hourly -> daily
timeinst_dt = datetime(time_inst','ConvertFrom','datenum'); 
flwinst_hr = timetable(timeinst_dt,along_inst');
flwinst_dy = retime(flwinst_hr, 'daily', @(x) mean(x, 'omitnan'));

acrsinst_hr = timetable(timeinst_dt,acrs_inst');
acrsinst_dy = retime(acrsinst_hr, 'daily', @(x) mean(x, 'omitnan'));


% Next, make timetable for cumulative flow
timecumu_dt = datetime(time_cumu,'ConvertFrom','datenum');
flwcumu_dy = timetable(timecumu_dt',flow_cumu'); 

% Then create overall daily timeseries - to slot both flow series into
timeall = (datetime(flwinst_dy.timeinst_dt(1).Year,1,1):caldays(1):datetime(flwinst_dy.timeinst_dt(end).Year,12,31)); % Master time array (daily, yearlong)

% Create empty NaN vectors for each flow type, then fill with the indices
% of the 'intersected' values with the master timeseries
install = NaN(1,length(timeall));
[~,iid] = intersect(timeall,flwinst_dy.timeinst_dt);
install(iid) = flwinst_dy.Var1;

cumuall = NaN(1,length(timeall));
[~,cid] = intersect(timeall,flwcumu_dy.Time);
cumuall(cid) = flwcumu_dy.Var1;

% Combine 'timeall' with the 'instflw' and 'cumuflw' fill-ins
flwsall_df = table(timeall',install',cumuall');
flwsall_df.Properties.VariableNames = {'Date' 'Inst_flow' 'Cumu_flow'};

writetable(flwsall_df,'NH10_Flows_Inst_Cumu.csv');

% %%% ALSO write table of Along and Across flows -> for GAMs
inst_df = table(flwinst_dy.timeinst_dt,flwinst_dy.Var1,acrsinst_dy.Var1);
inst_df.Properties.VariableNames{1} = 'Date';
inst_df.Properties.VariableNames{2} = 'Along_flow_inst';
inst_df.Properties.VariableNames{3} = 'Across_flow_inst';

writetable(inst_df,'NH10_AlongAcrossFlows_daily.csv');



