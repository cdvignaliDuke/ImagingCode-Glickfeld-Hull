function [output] = FA_Trl(Early,tCycleNum,trialnum,Leverup,Leverdown)
idx_early = find(Early==1);
output.FA_Num = 0; % release after the onset of cycle 3
output.fidget_Num = 0; % release before the onset of cycle 3

for i = 1: length(idx_early)
    idx_temp = [];
    cycle = [];
    trialoffs= [];
    Triallength = [];
    FA_RT = [];
    idx_temp = idx_early(i);% early trials
    cycle = tCycleNum(idx_temp);
    Triallength =Leverup(idx_temp)-Leverdown(idx_temp);    
    if cycle>=3  % release after cycle 3
        output.FA_Num =  output.FA_Num + 1;
    else
        output.fidget_Num = output.fidget_Num + 1;
    end
end
% FA rate is calculated here as percentage of total trial numbers


[c d]= binofit(output.FA_Num,trialnum);
output.FA_confi(1,1:2) = d;
output.FA(1,1)=c;


[c d]= binofit(output.fidget_Num,trialnum);
output.fidget_confi(1,1:2) = d;
output.fidget(1,1)=c;



end