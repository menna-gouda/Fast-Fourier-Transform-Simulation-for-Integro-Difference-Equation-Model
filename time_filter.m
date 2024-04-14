function [Isegs Tsegs]=time_filter(Isegs,Tsegs,time_threshold)

[p,q]=size(Isegs);
flag=[];
for i=1:p
   if Tsegs(i)<time_threshold
       flag=[flag 1];
   else 
       flag=[flag 0];
   end
end
 to_del=find(flag==1);
 Isegs(to_del,:)=[];
 Tsegs(to_del,:)=[];
end

