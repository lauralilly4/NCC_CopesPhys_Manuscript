%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         NOAA/OSU Project #1:        %%%
%%%  NH10 Alongshore Flow - Cumul Inds  %%%
%%%     S2: Calculate Cumul. Mags.    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From Brandy Cervantes (M. Fewings group)




%% ======= Step 2: Get yearly cumul. ALF =======
% Calculate *daily-avg* flows from all hourly values
hrflw = timetable(flwdt',flwalngavg', 'VariableNames', {'Alongshore_flow'}); % Get timetable of dates
% dyflwavg = groupsummary(hrflw, 'Time', 'day', 'mean')); % Calculate daily avgs
dyflw = retime(hrflw,'daily','mean'); % Calculate daily avgs
dyflwdt = dyflw.Time;
dyflwavg = dyflw.Alongshore_flow;


% Start for-loop
flwyrsunq = unique(dyflwdt.Year); % Get unique years from all dates

flwyrssub = flwyrsunq(ismember(flwyrsunq,wndyr)); % Get subset of ALF years that are also in ST/FT years (I don't love this code setup)
sprtrsub = wndst(ismember(wndyr,flwyrsunq),:); % Get subsets of Spring and Fall Trans years that match ALF years
faltrsub = wndft(ismember(wndyr,flwyrsunq),:);

yrflwmags = NaN(length(flwyrssub),365);

for yd = 2:length(flwyrssub) % Start w/ Year 2 (1998) because Year 1 (1997) isn't full
    yid = find(dyflwdt.Year == yrsunq(yd));
    % dchnk = dyflwdt(yid);
    fchnk = dyflwavg(yid);
    fsub = fchnk(sprtrsub(yd):faltrsub(yd)-1);
    fcumu = cumsum(fsub,'omitnan');
    yrmag = NaN(1,365);
    yrmag(sprtrsub(yd):faltrsub(yd)-1) = fcumu;
    yrflowmag(yd,:) = yrmag;
end

xaxdt = datetime(dyflwdt(yid),'Format','MMM-dd');
xaxdt(60) = []; % Remove Feb-29 (because this dt vector comes from a Leap Year)


%% ======= Step 3: Plot Cumul mags (alongshore flow) =======
figure(4);

hold on
for p=1:length(yrflowmag(:,1))
    plot(xaxdt(60:end),yrflowmag(p,60:end));
    text(xaxdt(302),min(yrflowmag(p,:)),num2str(flwyrssub(p)));
end
hold off
