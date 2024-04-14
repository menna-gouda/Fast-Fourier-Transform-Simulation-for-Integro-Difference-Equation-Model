%% Loading data 
data=["new_data.mat","meta_data.mat","IDindx.mat",'final_long.mat',...
    "m","f","a8","a20","fl","s","all_BICs","all_Parms","all_Flags","CIs"];
for i=1:length(data)
    load(data(i));
end
%% Estimate saturation rate "ab" and mean velocity "c"
% ab: 1/mean(travelling time)
ab=1/mean(final_long(:,2));  % 0.0590
% c: mean(travelling distance/travelling time)
c=mean(final_long(:,1)./final_long(:,2)); % 784.3065
%% Fitting 
r=final_long(:,1);  
alpha=mean(r); N=length(r); d=1;

format short g;
% nll modified bessel function of the second kind 
[parms1, fval1, iflag1]= fminsearch(@(x) nll_bessel(r, x), [alpha]) 
% nll type II
[parms2, fval2, iflag2]= fminsearch(@(x) nll_typeII(r, x), [alpha d]) 
% nll type III
[parms3, fval3, iflag3]= fminsearch(@(x) nll_typeIII(r, x), [alpha d]) 
% nll model 4
[parms4, fval4, iflag4]= fminsearch(@(x) nll_model4(r, x), [alpha d]) 
% nll model 4 new
[parms5, fval5, iflag5]= fminsearch(@(x) nll_model4_new(r, x), [alpha d]) 
%% Visualization
r=final_long(:,1); N=length(r);
[counts,centers]=hist(r,30); 
delta_r=(centers(2)-centers(1));  % bin width  
hist(r/1000,30)
format short g;

rc=linspace(1,71001,71001);
drc=rc(2)-rc(1);

% bessel NLL
alpha=parms1(1);
[p1,c1]=bessel(alpha);
hold on, plot(rc/1000,N*delta_r*p1*c1,"r"), 

% grad descent type II NLL
alpha=parms2(1);d=parms2(2);
[p2,c2]=typeII(parms2);
hold on, plot(rc/1000,N*delta_r*c2*p2,"m"), 

% grad descent type III NLL
alpha=parms3(1);d=parms3(2);
[p3,c3]=typeIII(parms3);
hold on, plot(rc/1000,N*delta_r*c3*p3,"g"), 

% grad descent model 4 NLL
alpha=parms4(1);d=parms4(2);
[p4,c4]=model4(parms4);
hold on, plot(rc/1000,N*delta_r*c4*p4,"c"), hold off

xlabel('Distance (km)'), ylabel('Frequency')
legend("Raw data", "K_1","K_2", "K_3", "K_4")
%% BIC Calculations
N=length(r);
format short g;
% bessel 
BIC1= 2*fval1 + 1*log(.5*N/pi)   % 4807.8
% grad type II
BIC2=  2*fval2 + 2*log(.5*N/pi)  % 4805.4
% grad type III
BIC3=  2*fval3 + 2*log(.5*N/pi)  % 4803.9
% grad model 4
BIC4=  2*fval4 + 2*log(.5*N/pi)  % 4828.1
% grad model 4 new
BIC5=  2*fval5 + 2*log(.5*N/pi)  % 4825
% Comparisons  
BIC13=exp(0.5*(BIC1-BIC3)); % 6.9748
BIC23=exp(0.5*(BIC2-BIC3)); % 2.0954
BIC43=exp(0.5*(BIC4-BIC3)); % 1.8349e+05
BIC53=exp(0.5*(BIC5-BIC3)); % 38525
%% Fitting partitioned data 
%JP:  note it's not just best but significantly best that matters
BIC_m=fit(m)     
[indx_best,best,indx_secbest,secbest]=bestBICs(BIC_m)   % best: 3, 2nd best: 2
BIC_f=fit(f)    
[indx_best,best,indx_secbest,secbest]=bestBICs(BIC_f)   % best: 1, 2nd best: 3 (Interesting!)
BIC_a8=fit(a8);      
[indx_best,best,indx_secbest,secbest]=bestBICs(BIC_a8)  % best: 3, 2nd best: 2
BIC_a20=fit(a20);    
[indx_best,best,indx_secbest,secbest]=bestBICs(BIC_a20) % best: 1, 2nd best: 3 (Interesting!)
%%
% these need to be changed to the new scheme (episode and yr)
BIC_s=fit(s);    
[indx_best,best,indx_secbest,secbest]=bestBICs(BIC_s)   % best: 2, 2nd best: 3
%%
BIC_fl=fit(fl);      
[indx_best,best,indx_secbest,secbest]=bestBICs(BIC_fl)  % best: 1, 2nd best: 3
%% Descriptive Statistics
% JP: some deer are multiple dispersers
unique_long=unique(final_long(:,3)); % 89 deer doing long dispersals
frac_long_deer=length(unique_long)/363;  % 0.2452
%%
all_males=find(meta_data{1:363,"sex"}=="Male");  % 185
frac_long_males=length(m)/length(all_males);     % 0.49189
all_females=find(meta_data{1:363,"sex"}=="Female");  % 178
frac_long_females=length(f)/length(all_females); % 0.1236
%% 
dead=find(meta_data{1:363,"reasonoff1"}=='mortality'); % 174 out of 363 died
frac_dead=length(dead)/363;  % 0.4793  Sad :(
%%
dead_long_dists=find(final_long(:,9)==0); % 48 out of the 113 dispersals the deer dies in the end
frac_dead_dists=length(dead_long_dists)/113; % 0.4248 (percentage of dispersals the deer dies in the end)
%%
dead_long_deer=unique(final_long(dead_long_dists,8)); % 43 long-dispesing deer died
%%
frac_male_deer=length(all_males)/363;              % 0.50964
frac_female_deer=length(all_females)/363;          % 0.49036
%% 
alive=363-dead; % the next part is repeated many times, it can be a function in the future
dead_males=find(meta_data{1:363,"sex"}=="Male" & meta_data{1:363,"reasonoff1"}=="mortality"); % 97 dead males 
alive_males=length(all_males)-length(dead_males); % 88
frac_males_alive=alive_males/length(all_males); % 0.4757
%%
dead_females=find(meta_data{1:363,"sex"}=="Female" & meta_data{1:363,"reasonoff1"}=="mortality"); % 77 dead females 
alive_females=length(all_females)-length(dead_females); % 101
frac_females_alive=alive_females/length(all_females);  % 0.5674
%%
all_age8=find(meta_data{1:363,"ageatcol1"}==8);    % 259
frac_long_a8=length(a8)/length(all_age8);          % 0.39768
all_age20=find(meta_data{1:363,"ageatcol1"}==20);  % 104
frac_long_a20=length(a20)/length(all_age20);       % 0.096154
%%
frac_a8_deer=length(all_age8)/363;                 % 0.7135
frac_a20_deer=length(all_age20)/363;               % 0.2865
%%
dead_a8=find(meta_data{1:363,"ageatcol1"}==8 & meta_data{1:363,"reasonoff1"}=="mortality"); % 122 dead 8mo
alive_a8=length(all_age8)-length(dead_a8);         % 137
frac_a8_alive=alive_a8/length(all_age8);           % 0.5290
%%
dead_a20=find(meta_data{1:363,"ageatcol1"}==20 & meta_data{1:363,"reasonoff1"}=="mortality"); % 52 dead 20mo
alive_a20=length(all_age20)-length(dead_a20);      % 52 
frac_a20_alive=alive_a20/length(all_age20);        % 0.5000
%%
mode(final_long(:,6)); % 12 
hist(final_long(:,6),20) % no long dispersals in the fifth year max(new_data(:,3)) = 1.5911e+03 (5th yr)
[n,bin] = hist(final_long(:,6),unique(final_long(:,6)));
[~,idx] = sort(-n);
%% 
n(idx)  % frequencies in descending order  24    19    16    13    13    12     8     8
transpose(bin(idx)) % corresponding values 21    12    32    11    31    42    41    22
%% Calculate # of deer dispersing at their first spring s1, first fall f1, or second spring s2
% I didn't use the mortality day as the collar could be lost or follow-up
% has stopped in any other way. So I used the dates the deer were tracked at only (from T in new_data).
% alive_at_season=zeros(10,1);
% code=[11 1; 21 2; 31 3; 41 4; 51 5; 12 6; 22 7; 32 8; 42 9; 52 10];
disp_season=zeros(363,3);
for n=1:363
    start_dy=new_data(IDindx(n,1),3);  
    end_dy=new_data(IDindx(n,2),3);    
    age=meta_data{n,"ageatcol1"};
    start_end_season(n,1)=calc_season(start_dy);
    start_end_season(n,2)=calc_season(end_dy);
    seasons_count=start_end_season(n,2)-start_end_season(n,1);
    collar_period=end_dy-start_dy;
    if collar_period>180  % a season is at least 180 days
        if seasons_count==0
            % seasons are s1, f1, and s2
            seasons(n,1)=1;
            seasons(n,2)=1;
            seasons(n,3)=1;
        elseif seasons_count==1 
            if age==8
                seasons(n,1)=1;
                seasons(n,2)=1;
            else
                seasons(n,3)=1;
            end
        elseif seasons_count==-1
            seasons(n,2)=1;
            seasons(n,3)=1;
        end
    else
        if age==8
            seasons(n,1)=1;
        else
            seasons(n,3)=1;
        end
    end
    % season_indx1=find(code(:,1)==start_end_season(n,1));
    % season_indx2=find(code(:,1)==start_end_season(n,2));
    % for i=season_indx1
    %     alive_at_season(i)=alive_at_season(i)+1;
    % end
    
    ind=find(final_long(:,8)==n);
    if length(ind)>1
        for i=ind(1):ind(length(ind))
            dispersal_start_day=final_long(i,7);
            period_till_dispersal=dispersal_start_day-start_dy;
            dispersal_season=calc_season(dispersal_start_day);
            if dispersal_season==2
                disp_season(n,2)=disp_season(n,2)+1;
           elseif dispersal_season==1 && period_till_dispersal<180
                if age==8
                    disp_season(n,1)=disp_season(n,1)+1;
                else
                    disp_season(n,3)=disp_season(n,3)+1;
                end
            elseif dispersal_season==1 && period_till_dispersal>180
                disp_season(n,3)=disp_season(n,3)+1;
            end
        end 
    elseif length(ind)==1
        dispersal_start_day=final_long(ind,7);
        period_till_dispersal=dispersal_start_day-start_dy;
        dispersal_season=calc_season(dispersal_start_day);
        if dispersal_season==2
            disp_season(n,2)=1;
       elseif dispersal_season==1 && period_till_dispersal<180
            if age==8
                disp_season(n,1)=1;
            else
                disp_season(n,3)=1;
            end
        elseif dispersal_season==1 && period_till_dispersal>180
            disp_season(n,3)=1;
        end
    end
end  % This way we are not capturing if somedeer dispersed four times or twice a season.  
%% Bar graph showing dispersals to total alive population per season
format short g;
s1=sum(disp_season(:,1))*100/sum(seasons(:,1));  % 13.31  % 13.312
f1=sum(disp_season(:,2))*100/sum(seasons(:,2));  % 16.72  % 18.033
s2=sum(disp_season(:,3))*100/sum(seasons(:,3));  %  5.68  % 7.4236
x=categorical({'First Spring', 'First Fall', 'Second Spring'});
x=reordercats(x,{'First Spring', 'First Fall', 'Second Spring'});
b=bar(x,[s1;f1;s2]);
b.FaceColor='flat';
text(1,s1,num2str(s1'),'vert','bottom','horiz','center'); 
text(2,f1,num2str(f1'),'vert','bottom','horiz','center'); 
text(3,s2,num2str(s2'),'vert','bottom','horiz','center'); 
box off
b.CData(3,:)=[.2 .6 .5];
b.CData(2,:)=[.5 0 .5];
ylabel("Percentage of long-distance dispersals to total alive deer")
title("Long-distance dispersals at each season")
%% kinda old stuff here
freq_long_season=[11 13; 21  24; 31  13; 41  8; 51 0; 12 19; 22 8; 32 16; 42 12; 52 0];
frac_long_alive=freq_long_season(:,2)./alive_at_season;
% 0.158 0.12 0.04 0.024 0 0.10 0.04 0.13 0.18 NaN
%   1     2   3     4   5   6    7   8    9   10
% hist(frac_long_alive)
plot(sort(freq_long_season(:,1)),frac_long_alive)
xlabel("year,season"); ylabel("fraction of alive deer doing long dispersals");
%% old as well
bar_data=[13 19; 24 8; 13 16; 8 12; 0 0];
bar_fracs=[frac_long_alive(1) frac_long_alive(6); frac_long_alive(2) frac_long_alive(7); ...
    frac_long_alive(3) frac_long_alive(8); frac_long_alive(4) frac_long_alive(9); ...
    frac_long_alive(5) frac_long_alive(10)];
bar_x=2017:1:2021;
figure, bar(bar_x,bar_data)
legend("episode 1","episode 2")
title("count of dispersers at each season")
figure, bar(bar_x,bar_fracs)
title("fraction of alive dispersers at each season")
legend("episode 1","episode 2")
%% Bootstrapping
all_BICs=[];all_Parms=[];all_Flags=[];
while sum(all_Flags)~=1000 
    rand_idx=ceil(113*rand(113,1));
    [BICs,Parms,Flags]=fit(final_long(rand_idx,1)); 
    if (Flags==1) % only if all four == 1
       all_BICs=[all_BICs;BICs]; 
       all_Parms=[all_Parms;Parms];
       all_Flags=[all_Flags;Flags];
    end
end
%% Boostrap histograms
nominal_BICs=[BIC1,BIC2,BIC3,BIC4];
nominal_Parms=[parms1,parms2,parms3,parms4];
n=1;
for i=1:4
    create_bootstrap_hists(all_BICs(:,i),nominal_BICs(i),i,"BIC");
    create_bootstrap_hists(all_Parms(:,n),nominal_Parms(n),i,text('string',"$\alpha$",'interpreter','latex'));
    n=n+1;
    if i > 1
        create_bootstrap_hists(all_Parms(:,n)+1,nominal_Parms(n)+1,i,text('string',"$\beta$",'interpreter','latex'));
        n=n+1;
    end
end
%% fitting males and females to get nominal values for them
% This could have been more efficient by using the function fit()
% using pre-made m, f vectors
alpha_m=mean(m); N_m=length(m); d_m=1;
alpha_f=mean(f); N_f=length(f); d_f=1;

format short g;
% nll modified bessel function of the second kind 
[parms1_m, fval1_m, iflag1_m]= fminsearch(@(x) nll_bessel(m, x), [alpha_m]) 
[parms1_f, fval1_f, iflag1_f]= fminsearch(@(x) nll_bessel(f, x), [alpha_f]) 
% nll type II
[parms2_m, fval2_m, iflag2_m]= fminsearch(@(x) nll_typeII(m, x), [alpha_m d_m]) 
[parms2_f, fval2_f, iflag2_f]= fminsearch(@(x) nll_typeII(f, x), [alpha_f d_f]) 
% nll type III
[parms3_m, fval3_m, iflag3_m]= fminsearch(@(x) nll_typeIII(m, x), [alpha_m d_m]) 
[parms3_f, fval3_f, iflag3_f]= fminsearch(@(x) nll_typeIII(f, x), [alpha_f d_f]) 
% nll model 4
[parms4_m, fval4_m, iflag4_m]= fminsearch(@(x) nll_model4(m, x), [alpha_m d_m]) 
[parms4_f, fval4_f, iflag4_f]= fminsearch(@(x) nll_model4(f, x), [alpha_f d_f]) 

% bessel 
BIC1_m= 2*fval1_m + 1*log(.5*N_m/pi) % NLL
BIC1_f= 2*fval1_f + 1*log(.5*N_f/pi) % NLL
% grad type II
BIC2_m=  2*fval2_m + 2*log(.5*N_m/pi)  % NLL
BIC2_f=  2*fval2_f + 2*log(.5*N_f/pi)  % NLL
% grad type III
BIC3_m=  2*fval3_m + 2*log(.5*N_m/pi)  % NLL
BIC3_f=  2*fval3_f + 2*log(.5*N_f/pi)  % NLL
% grad model 4
BIC4_m=  2*fval4_m + 2*log(.5*N_m/pi)  % NLL
BIC4_f=  2*fval4_f + 2*log(.5*N_f/pi)  % NLL
%% Bootstrapping again (while accounting for differences in sex)
all_BICs_j=[];all_Parms_j=[];all_Flags_j=[];  % joint sample
all_BICs_m=[];all_Parms_m=[];all_Flags_m=[];  % male sample
all_BICs_f=[];all_Parms_f=[];all_Flags_f=[];  % female sample

while (sum(all_Flags_j)~=1000) | (sum(all_Flags_m)~=1000) | (sum(all_Flags_f)~=1000)
    % disp("hello");
  
    rand_idx=ceil(113*rand(113,1));
    final_dists=final_long(rand_idx,:);

    % males
    male_indicies=find(final_dists(:,4)==0); 
    male_dists=final_dists(male_indicies,:);        
    m_btstrp=male_dists(:,1);
    % females
    female_indicies=find(final_dists(:,4)==1); 
    female_dists=final_dists(female_indicies,:);   
    f_btstrp=female_dists(:,1);

    [BICs_j,Parms_j,Flags_j]=fit(final_dists(:,1));
    [BICs_m,Parms_m,Flags_m]=fit(m_btstrp); 
    [BICs_f,Parms_f,Flags_f]=fit(f_btstrp);
    
    % ignoring samples if their fit doesn't converge or BIC goes complex
    if (Flags_j==1) & (Flags_m==1) & (Flags_f==1) & (isreal(BICs_j)) & (isreal(BICs_m)) & (isreal(BICs_f))     
       % joint sample
       all_BICs_j=[all_BICs_j;BICs_j]; 
       all_Parms_j=[all_Parms_j;Parms_j];
       all_Flags_j=[all_Flags_j;Flags_j];
       % male sample
       all_BICs_m=[all_BICs_m;BICs_m]; 
       all_Parms_m=[all_Parms_m;Parms_m];
       all_Flags_m=[all_Flags_m;Flags_m];
       % female sample
       all_BICs_f=[all_BICs_f;BICs_f]; 
       all_Parms_f=[all_Parms_f;Parms_f];
       all_Flags_f=[all_Flags_f;Flags_f];
    end
end
%% Bootstrap histograms again (for males and females version)
nominal_BICs_j=[BIC1,BIC2,BIC3,BIC4];
nominal_Parms_j=[parms1,parms2,parms3,parms4];
nominal_BICs_m=[BIC1_m,BIC2_m,BIC3_m,BIC4_m];
nominal_Parms_m=[parms1_m,parms2_m,parms3_m,parms4_m];
nominal_BICs_f=[BIC1_f,BIC2_f,BIC3_f,BIC4_f];
nominal_Parms_f=[parms1_f,parms2_f,parms3_f,parms4_f];
CIs_BICs_j=CIs(1:4,:); CIs_BICs_m=CIs(5:8,:); CIs_BICs_f=CIs(9:12,:); 
CIs_Parms_j=CIs(13:19,:); CIs_Parms_m=CIs(20:26,:); CIs_Parms_f=CIs(27:33,:);
CIs_delta_BICs_j=CIs(34:37,:); CIs_delta_BICs_m=CIs(38:41,:); CIs_delta_BICs_f=CIs(42:45,:); 
n=1;
for i=1:4
    %create_bootstrap_hists(all_BICs_j(:,i),nominal_BICs_j(i),i," BIC for joint sample using model ",CIs_BICs_j(i,1),CIs_BICs_j(i,2));
    %create_bootstrap_hists(all_BICs_m(:,i),nominal_BICs_m(i),i," BIC for male sample using model ",CIs_BICs_m(i,1),CIs_BICs_m(i,2));
    %create_bootstrap_hists(all_BICs_f(:,i),nominal_BICs_f(i),i," BIC for female sample using model ",CIs_BICs_f(i,1),CIs_BICs_f(i,2));
    %create_bootstrap_hists(all_Parms_j(:,n),nominal_Parms_j(n),i," Parameter $\alpha$ for joint sample model using model ",CIs_Parms_j(n,1),CIs_Parms_j(n,2));
    %create_bootstrap_hists(all_Parms_m(:,n),nominal_Parms_m(n),i," Parameter $\alpha$ for male sample model using model ",CIs_Parms_m(n,1),CIs_Parms_m(n,2));
    %create_bootstrap_hists(all_Parms_f(:,n),nominal_Parms_f(n),i," Parameter $\alpha$ for female sample model using model ",CIs_Parms_f(n,1),CIs_Parms_f(n,2));
    create_bootstrap_hists(all_BICs_j(:,i)-all_BICs_j(:,3),nominal_BICs_j(i)-nominal_BICs_j(3),i," $\Delta$ BIC for joint sample for model ",CIs_delta_BICs_j(i,1),CIs_delta_BICs_j(i,2));
    %create_bootstrap_hists(all_BICs_m(:,i)-all_BICs_m(:,2),nominal_BICs_m(i)-nominal_BICs_m(2),i," $\Delta$ BIC for male sample for model ",CIs_delta_BICs_m(i,1),CIs_delta_BICs_m(i,2));
    %create_bootstrap_hists(all_BICs_f(:,i)-all_BICs_f(:,1),nominal_BICs_f(i)-nominal_BICs_f(1),i," $\Delta$ BIC for female sample for model ",CIs_delta_BICs_f(i,1),CIs_delta_BICs_f(i,2));
    n=n+1;
    if i > 1
       %create_bootstrap_hists(all_Parms_j(:,n)+1,nominal_Parms_j(n)+1,i," Parameter $\beta$ for joint sample model ",CIs_Parms_j(n,1),CIs_Parms_j(n,2));
       %create_bootstrap_hists(all_Parms_m(:,n)+1,nominal_Parms_m(n)+1,i," Parameter $\beta$ for male sample model ",CIs_Parms_m(n,1),CIs_Parms_m(n,2));
       %create_bootstrap_hists(all_Parms_f(:,n)+1,nominal_Parms_f(n)+1,i," Parameter $\beta$ for female sample model ",CIs_Parms_f(n,1),CIs_Parms_f(n,2));
       n=n+1;
    end
end
%% Deciding which is the best model for each group: joint, male, and female
best_model=[0 0 0];
for j=1:3
    if j==1
        myBIC=all_BICs_j;
    elseif j==2
        myBIC=all_BICs_m;
    else 
        myBIC=all_BICs_f;
    end
    min_sum=inf;
    for i=1:4
        if mean(myBIC(:,i))<min_sum  % should I use sum?! or mode? or mean?
           min_sum=mean(myBIC(:,i));
           best_model(j)=i;   % [3 2 1] using sum 
        end                   % [3 2 3] using mode
    end                       % [3 2 1] using mean
end
%% Comparing BICs male + female vs joint
% map = brewermap(2,'Set1'); 
combined_BICs=all_BICs_m(:,2)+all_BICs_f(:,1);
figure, 
h1=histogram(all_BICs_j(:,3),'facecolor',[0 0 1],'edgecolor','none'), hold on
h2=histogram(combined_BICs,'facecolor',[1 0 0],'edgecolor','none'), hold off
title('Distribution of BIC values for male-female and joint samples','interpeter','latex')
xl=xline(BIC3, '--k',strcat("Nominal value: ",num2str(round(BIC3))));
xll=xline(BIC2_m+BIC1_f, '--k',strcat("Nominal value: ",num2str(round(BIC2_m+BIC1_f))));
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'right';
xl.Color=[0 0 1];
xl.LineWidth=1.3;
xll.Color=	[1 1 0];
xll.LineWidth=1.3;
xll.LabelVerticalAlignment = 'bottom';
xll.LabelHorizontalAlignment = 'left';
xlabel('BIC values','interpeter','latex'); ylabel('frequency','interpeter','latex')
legend("Joint","Male-Female","joint", "Male-Female")
%% getting delta BICs for male+female vs joint
combined_BICs=all_BICs_m(:,2)+all_BICs_f(:,1);
figure, 
h1=histogram(all_BICs_j(:,3)-combined_BICs,'facecolor',[0 0 1],'edgecolor','none')
title('Distribution of \Delta BIC for male-female vs. joint sample')
xl=xline(BIC3-BIC2_m-BIC1_f, '--k',strcat("nominal value: ",num2str(BIC3-BIC2_m-BIC1_f)));
% xll=xline(BIC2_m+BIC1_f, '--k',strcat("male-female sample nominal value: ",num2str(round(BIC2_m+BIC1_f))));
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'right';
xl.LineWidth=1.3;
% xl.Color=[0 1 1];
xlabel('\Delta BIC values'); ylabel('frequency')
%%
combined_BICs=all_BICs_m(:,2)+all_BICs_f(:,1);
delta_BICs_mfj=all_BICs_j(:,3)-combined_BICs;
ind_delta_BIC_mfj=find(delta_BICs_mfj>0);
length(ind_delta_BICs_mfj)/10
%% Plotting settling rate functions
a=5;b=3; t=linspace(0,40);
h1=a*b; h2=(a*b.*t)./(a+t); h3=(a*b.*t.^2)./(a^2+t.^2); h4=(b.*t.^2)./(a+t);
plot(t,h1+t*0,t,h2,t,h3,t,h4,'LineWidth',1.3)
ylim([0 40]);
yline(15,'-','h(t) = ab')
xline(5,'--','t = a')
legend("h_1","h_2","h_3","h_4")
xlabel('time (days)')
ylabel('Settling rate h(t)')
title("Settling rate functions")
% text(38,16,"ab")
% text(5.2,1,"a")
%% Plotting settling rate functions in another way
a = 5;
b = 3;
t = linspace(0, 40);

h1 = a * b;
h2 = (a * b .* t) ./ (a + t);
h3 = (a * b .* t.^2) ./ (a^2 + t.^2);
h4 = (b .* t.^2) ./ (a + t);

figure;
plot(t, h1 + t * 0, 'k', 'LineWidth', 1.3); % Plotting h1 in black
hold on;
plot(t, h2, 'k--', 'LineWidth', 1.3); % Plotting h2 in black dashed line
plot(t, h3, 'k-.', 'LineWidth', 1.3); % Plotting h3 in black dash-dot line
plot(t, h4, 'k:', 'LineWidth', 1.3); % Plotting h4 in black dotted line
ylim([0 40]);
yline(15, 'k', 'h(t) = ab'); % Adding a horizontal line at y=15 with black color
xline(5, 'k--', 't = a'); % Adding a vertical line at x=5 with black dashed line
xlabel('Time (days)');
ylabel('Settling rate h(t)');
title('Settling rate functions');

% Adding legend with specified labels
legend('h_1', 'h_2', 'h_3', 'h_4', 'Location', 'NorthEast');

% Alternatively, you can add legend labels with LaTeX formatting
% legend('h_1', 'h_2', 'h_3', 'h_4', 'Interpreter', 'latex', 'Location', 'NorthEast');
%% Plotting trajectories, home and travel ranges
% figure(1),plot_traj(n) % whole trajectory of nth deer
%  for j=1:length(Tsegs)
%      figure(j+1),plot_home(n,j,Isegs) % jth home range of the nth deer
%  end
%  for k=1:length(dsegs)
%      figure(k+j+1),plot_trvl(n,k,dsegs) % kth travelling range of the nth deer
%  end
%% Draft plots for histogram distribtion and bessel function
% hist(r,nbins)
% R=[centers-delta_r/2 max(centers)+delta_r/2]; % bin edges
% xlabel('Distance (m)'), ylabel('Frequency')
% K=besselk(0,R/alpha); % Does this now represent the distribution of rs?
% plot(R/alpha,K) % use R instead of r, bec r is not ordered  
% How I know this is a pdf? % And how after integration we have a CDF? 
%% Draft
% bootstrpd_data=[];
% random_indicies=round(1+112.*rand(113,1500));
% [rows,colmns]=size(random_indicies);
% for i=1:rows
%     for j=1:colmns
%         bootstrpd_data(i,j)=final_long(random_indicies(i,j),1);
%     end
% end
%     [BICs,Parms,Flags]=fit(bootstrpd_data(:,k));
%% Deciding how to separate empirical fall and spring episodes
[N,edges]=histcounts(mod(final_long(:,7),365),30);
sum(N(1:15));
sum(N(17:30));
hist(mod(final_long(:,7),365),30) % 180 days is a good number (175-187) (172.5-184)
xlabel("time (days)"); ylabel("frequency of long-distance dispersal")
title("Episodes of juvenile long-distance dispersal")
%% Confidence Intervals
% all_BICs_j,all_BICs_m,all_BICs_f,all_Parms_j,all_Parms_m,all_Parms_f;
% CIs_BICs_j=CIs(1:4,:); CIs_BICs_m=CIs(5:8,:); CIs_BICs_f=CIs(9:12,:); 
% CIs_Parms_j=CIs(13:19,:); CIs_Parms_m=CIs(20:26,:); CIs_Parms_f=CIs(27:33,:);
% CIs_delta_BICs_j=CIs(34:37,:) % relative to model III;
% CIs_delta_BICs_m=CIs(38:41,:) % relative to model II;
% CIs_delta_BICs_f=CIs(42:45,:) % relative to model I;
for i=3:2:7
    [CI_N,CI_x]=hist(all_Parms_f(:,i)+1,100000);
    freq=cumsum(CI_N)/sum(CI_N);
    low_val=max(freq(freq<=0.05));
    high_val=min(freq(freq>=0.95));
    lower_index=max(find(freq==low_val));
    higher_index=min(find(freq==high_val));
    CIs(i+26,1)=CI_x(lower_index);
    CIs(i+26,2)=CI_x(higher_index);
end
%% Average mean dispersal distance
format short g
av_male_dists=mean(m);
av_female_dists=mean(f);
av_fall_dists=mean(fl);
av_spring_dists=mean(s);
av_age_8=mean(a8);
av_age_20=mean(a20);
%% Draft
format short g
mean(all_Parms_j(:,6))
format long g
mean(all_BICs_m(:,3))
%%

