function season=calc_season(day)

% episode 1 is spring 
% epidsode 2 is fall
if(day<180)  % doing it the exhaustive way :/
    season=1;               % episode 1 year 1
elseif(day>365 && day<545)
    season=1;               % episode 1 year 2
elseif(day>730 && day<910)
    season=1;               % episode 1 year 3
elseif(day>1095 && day<1275)
    season=1;               % episode 1 year 4
elseif(day>1460 && day<1640)
    season=1;               % episode 1 year 5
elseif(day>180 && day<=365)
    season=2;               % episode 2 year 1
elseif(day>545 && day<=730)
    season=2;               % episode 2 year 2
elseif(day>910 && day<=1095)
    season=2;               % episode 2 year 3
elseif(day>1275 && day<=1460)
    season=2;               % episode 2 year 4
elseif(day>1640 && day<=1825)
    season=2;               % episode 2 year 5
end
