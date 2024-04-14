function mean_std=calc_std(new_data,IDindx,days2filter,speed_threshold)

% This function calculates the average std of home range
std_d=[];
for n=1:363 
    x=new_data(IDindx(n,1):IDindx(n,2),1);
    y=new_data(IDindx(n,1):IDindx(n,2),2);
    T=new_data(IDindx(n,1):IDindx(n,2),3);
    dx=new_data(IDindx(n,1):IDindx(n,2),5);
    dy=new_data(IDindx(n,1):IDindx(n,2),6);
    % Calling the fourier function
    [Isegs, Tsegs]=FT_traj(T,dx,dy,days2filter,speed_threshold); 
    % Find average stdiance across home range
    [p,q]=size(Isegs);
    for i=1:p
        if Tsegs(i)>14  % filtering non home ranges
            new_x=x(Isegs(i,1):Isegs(i,2));
            new_y=y(Isegs(i,1):Isegs(i,2));
            mean_x=mean(new_x);
            mean_y=mean(new_y);
            new_std=sqrt(sum((new_x-mean_x).^2+(new_y-mean_y).^2)./(length(new_x)-1)); % calculating std in 2D
            std_d=[std_d; new_std]; 
        end
    end
mean_std=mean(std_d,'omitnan'); % averaging all std
end