%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            NOAA/OSU Project #1:             %%%
%%%     NH10 Alongshore Flow - Cumul Inds       %%%
%%%   v2 - just using flow values calculated by BTC %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From Brandy Cervantes (M. Fewings group)



%% ======= Step 1: Load files and assign variables =======
wndflwfl = load('Phys_files_2024/Wind_Indices_Vel_45N_V2_1JanCumFlow.mat');


% %% Assign variables
% Flow
flwcumu = wndflwfl.cum_vel;
flwtm = wndflwfl.time_vel;
flwdt = datetime(flwtm,'ConvertFrom','datenum');
flwyrs = unique(flwdt.Year);
flwyrsall = flwyrs(1):1:flwyrs(end); % Create a sequence of ALL years encapsulated by the timeseries (because three years 
                                     % DNE in the timeseries above)

% Winds (NEWPORT)
wndyr = wndflwfl.years';
wndtm = wndflwfl.wind_nwpo3.time;
wnddt = datetime(wndtm,'ConvertFrom','datenum');
wndst = wndflwfl.wind_nwpo3.st';
wndft = wndflwfl.wind_nwpo3.ft';

wndyrtrns = table(wndyr',wndst',wndft');
wndyrtrns.Properties.VariableNames = {'Year' 'Spr_YrDy' 'Fal_YrDy'};

% %% Save Transition Dates (from Wind File) -> to .csv
% writetable('wndyrtrns','NH10_ALF_TrnsDts.csv');


% Convert Year-Days to actual dates (by year)
sprtr_dts = [];
faltr_dts = [];

for yd=30:length(wndst) % Start w/ 1996 (idx no. = 30)
    yrd = wndyr(yd);
    sprdt = datetime(['1-Jan-',num2str(yrd)])+wndst(yd)-1;
    faldt = datetime(['1-Jan-',num2str(yrd)])+wndft(yd)-1;
    
    sprtr_dts = [sprtr_dts;sprdt];
    faltr_dts = [faltr_dts;faldt];
end



%% ======= Step 2: Create and assign a daily 'master' timeseries =======
timechnks = NaN(length(flwyrsall),365);
yearchnks = NaN(length(flwyrsall),365);

for yd = 1:length(flwyrsall)
    tchnk = flwdt(find(flwdt.Year == flwyrsall(yd)));
    ychnk = flwcumu(find(flwdt.Year == flwyrsall(yd)));
        
    if length(tchnk) == 0 % Check if dates exist for a year -> if not, create a new date-chunk
        tchnk = datetime(flwyrsall(yd),1,1):caldays(1):datetime(flwyrsall(yd),12,31);
    end
    dnmchnk = datenum(tchnk); % Either way, convert 'datetime' chunk -> 'datenum'
    
    if length(ychnk) == 0
        ychnk = NaN(1,365);
    end
  
    tmstr = (datetime(flwyrsall(yd),1,1):caldays(1):datetime(flwyrsall(yd),12,31)); % Master time array (daily, yearlong)
    dnmmstr = datenum(tmstr);
    ymstr = NaN(1,length(tmstr));
    [~,sid] = intersect(tmstr,tchnk);
    dnmmstr(sid) = dnmchnk;
    ymstr(sid) = ychnk;
    
    % Remove Leap Days, if exist
    if length(dnmmstr) == 366
        dnmmstr(60) = [];
        ymstr(60) = [];
    end
    
    timechnks(yd,:) = dnmmstr;
    yearchnks(yd,:) = ymstr;
    
end


% Master timeseries -> for plots
tmstrplt = (datetime(flwyrs(2),1,1):caldays(1):datetime(flwyrs(2),12,31)); % Master time array (daily, yearlong)

% Convert 'timechnks' -> 'datetime' format
dtchnks = datetime(timechnks,'ConvertFrom','datenum','Format','dd-MMM-yyyy');



%% ======= Step 3: Plot Cumul mags (alongshore flow) =======
figure(4);

hold on
for p=1:length(yearchnks(:,1))
    plot(tmstrplt,yearchnks(p,:));
    % text(xaxdt(302),min(yrflowmag(p,:)),num2str(flwyrssub(p)));
end
hold off



%% ======= Step 4: Export yearly cumul mags as .csv =======
% Reshape flow and time 'chunks' to 1D arrays
dtrshp = reshape(dtchnks',[],1);
flwrshp = reshape(yearchnks',[],1);

flwtbl = table(dtrshp,flwrshp);

writetable(flwtbl,'NH10_ALF_DyMag.csv');





% %%% Addendum #1, cont: Get 'instflw_dy' (from 'Fig6_S1_ALF_DataIn.m'
% code) and 'flwcumu' on same timeframe so they can be saved in one .csv
% file -> to be worked with vs. copepod analyses
tmstr_cmb = (datetime(instflw_dy.flwtim_dt(1).Year,1,1):caldays(1):datetime(instflw_dy.flwtim_dt(end).Year,12,31)); % Master time array (daily, yearlong)

% Create empty NaN vectors for each flow type, then fill with the indices
% of the 'intersected' values with the master timeseries
instmstr = NaN(1,length(tmstr_cmb));
[~,iid] = intersect(tmstr_cmb,instflw_dy.flwtim_dt);
instmstr(iid) = instflw_dy.along_filt;

cumumstr = NaN(1,length(tmstr_cmb));
[~,cid] = intersect(tmstr_cmb,flwdt);
cumumstr(cid) = flwcumu;

% Combine 'timemaster' with the 'instflw' and 'cumuflw' fill-ins
flwsall_df = table(tmstr_cmb',instmstr',cumumstr');
flwsall_df.Properties.VariableNames = {'Date' 'Inst_flow' 'Cumu_flow'};

writetable(flwsall_df,'NH10_Flows_Inst_Cumu.csv');
