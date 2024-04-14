function [Isegs, Tsegs] = FT_traj(T, dx, dy, days2filter, speed_threshold)
%FT_TRAJ This is a function which uses Fourier filtering to find signficant 
%   excursions from home range based on GPS collar data (based on WTD
%   collar data from upper Midwest).  
%
%   Jim Powell, Dec. 23, 2021,  jim.powell@usu.edu
%
%   INPUTS (all vectors of the same length):
%       T           vector of location times, in Julian days
%       dx          W->E changes in location at each time, meters 
%       dy          S->N changes in location at each time, meters 
%
%   CONTROL PARAMETERS (scalars)
%       days2filter         time window of filter, days
%                               -- duration of a `week'
%       speed_threshold     average speed threshold for big move, meters/day
%       
%   OUTPUTS
%       
%       Isegs       two columns which are starting, ending indices for the
%                   portion of trajectories which stay below the filtering
%                   threshold (i.e. bracketing potential home ranges)
%       Tsegs       time individual spent in potential home range (useful
%                   for further filtering)
%

%   Calculate cumulative changes from initial location
cumx=cumsum(dx); cumy=cumsum(dy);
%   Total square displacement from original location
cumd=sqrt(cumx.^2+cumy.^2);

%   Find indices so that averaging doesn't fall out of data window
NT=length(T);
daze=days2filter/2;
idxwk_min=find(T-T(1)>daze, 1);
idxwk_max=find(T(NT)-T>daze, 1, 'last');
Ts=T(idxwk_min:idxwk_max);
N=idxwk_max-idxwk_min+1;    % length of smoothed trajectory, subtracting a 
                            % half week on each end

%   Sample velocities over a week (initial smoothing)
v_smooth=zeros(N,1);
for n=1:N
    % Find indices for beginning and end of weekly filter window (times
    % steps not equal!)
    iwk_left=find(T-T(n+idxwk_min-1) > -daze,1);  %find relative first index in half week interval to left
    iwk_right=min(NT-1,find(T-T(n+idxwk_min-1) > daze,1)); % find last in half week on right
    deltat=T(iwk_right)-T(iwk_left);           % actual time interval ~ 1 week
    v_smooth(n)=(cumd(iwk_right)-cumd(iwk_left))/deltat;
end

%   Set up for Fourier smoothing
%   First, figure out number of modes to keep based on filtering window and
%   calculate a mask of 0/1s to keep those modes
Nsmooth=int16(ceil(4*(Ts(N)-Ts(1))/days2filter)); % modes to keep; window duration/length of week
%       is freq. number for a wave one `week' long.  Addtl factor of two because 
%       of Nyquist sampling and we want to resolve bumps of window duration 
%       and another factor of two because we will reflect the series to minimize 
%       oscillations at beginning and end
 
if (rem(Nsmooth,1)==0 && Nsmooth > 0) 
    mask=zeros(2*N,1);
    mask(1:Nsmooth+1)=1; 
    mask(2*N-Nsmooth:2*N)=1;    % freqs to keep are in initial and final positions of FFT vector
    
    %   Use FFT to generate smoothed data in wave space
    fv=fft([v_smooth; v_smooth(N:-1:1)]);   % reflect data
    fv2=fv.*mask;           % filter out high frequencies
    v_s2=real(ifft(fv2));   % invert the FFT
    v_s=v_s2(1:N);          % smoothed velocities; keep first half of reflected results
    
    %   Find indices when smoothed velocites exceed threshold
    idx_normal=find(abs(v_s)<speed_threshold);   % indices corresponding to home range moves
    idx_jump=find(diff(idx_normal)>1);      % indices of idx_normal at which there 
    %       is a jump between ranges, i.e. indices when rapid average movement
    %       begins
         
    %   Set up output variables
    Nsegs=length(idx_jump)+1;   % number of distinct home range segments
    Tsegs=zeros(Nsegs,1);       % to contain the temporal duration of segments
    Isegs=zeros(Nsegs,2);       % start and end of relevant indices for each segment
    idx_hr=Isegs;               % working variable to keep track of home range indices 
                                % relative to original inputs
    
    if (Nsegs==1)       % no jumps, individual stays at home
        idx_hr(1,1)=0; idx_hr(1,2)=idxwk_max-1;
    elseif (Nsegs==2)   % one jump, before/after ranges
        idx_hr(1,1)=0; idx_hr(1,2)=idxwk_min+idx_normal(idx_jump(1));
        idx_hr(2,1)=idxwk_min+idx_normal(idx_jump(1)+1); idx_hr(2,2)=idxwk_max-1;
    else                % multiple jumps among ranges
        idx_hr(1,1)=0; idx_hr(1,2)=idxwk_min+idx_normal(idx_jump(1));
        idx_hr(Nsegs,1)=idxwk_min+idx_normal(idx_jump(Nsegs-1)+1); 
        idx_hr(Nsegs,2)=idxwk_max-1;
        for j=2:Nsegs-1
            idx_hr(j,1)=idxwk_min+idx_normal(idx_jump(j-1)+1);
            idx_hr(j,2)=idxwk_min+idx_normal(idx_jump(j));
        end
    end
    Isegs=idx_hr+1;    % starting index is 1, not zero
    
    % calculate the duration of residence in potential range segments
    for j=1:Nsegs
        Tsegs(j)=T(Isegs(j,2))-T(Isegs(j,1));
    end
else 
    disp("no");
    Isegs=[];
    Tsegs=[];
end

end

