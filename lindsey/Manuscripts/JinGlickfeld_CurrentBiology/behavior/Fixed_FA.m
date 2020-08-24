function [output] = Fixed_FA(Early,input,tCycleNum,trialnum,Leverup,Leverdown)

idx_early = find(Early==1);
output.FA_Num = 0;
output.S_FA_Num = 0;
output.L_FA_Num = 0;

output.FA_RT=[];
output.S_FA_RT = [];
output.L_FA_RT = [];
output.cycle = [];



for i = 1: length(idx_early)
    idx_temp = [];
    cycle = [];
   
    Triallength = [];
    FA_RT = [];
    idx_temp = idx_early(i);% early trials
    cycle = tCycleNum(idx_temp);
    Triallength =Leverup(idx_temp)-Leverdown(idx_temp); 
    FA_RT = (Triallength-250.*(cycle-1)-(double(input.stimOnTimeMs)*(cycle-1)));
    if cycle>=3  % only consider FA occur after cycle 2, because target never occur on cycle 1&2; through out early release for target
        %calculate collapsed FAs
        
            % if release in the RT window 200-550 ms 
            if  FA_RT>=200 && FA_RT<=550
              
                output.FA_Num = output.FA_Num+1;
                output.FA_RT = [output.FA_RT, FA_RT];
                output.cycle = [output.cycle, cycle];
                
                if cycle<=5   % cycle number 3,4,5 as short cycles
                    output.S_FA_Num = output.S_FA_Num+1;
                    output.S_FA_RT = [output.S_FA_RT, FA_RT];
                end
                if cycle>=7 && cycle<=9  % cycle number 7,8,9 as long cycles
                    output.L_FA_Num = output.L_FA_Num+1;
                    output.L_FA_RT = [output.L_FA_RT, FA_RT];
                end
                
            end
            % if release very fast and previous off time is 250ms off, then
            % this is count as FA for the previous cycle..
            
            if FA_RT<200 && cycle>3
                output.FA_Num = output.FA_Num+1;
                output.FA_RT = [output.FA_RT, (FA_RT+250+double(input.stimOnTimeMs))];
                output.cycle = [output.cycle, cycle-1];
                
                if cycle<=6   % cycle number 3,4,5 as short cycles
                    output.S_FA_Num = output.S_FA_Num+1;
                    output.S_FA_RT = [output.S_FA_RT, (FA_RT+250+double(input.stimOnTimeMs))];
                end
                if cycle>=8 && cycle<=10  % cycle number 7,8,9 as long cycles
                    output.L_FA_Num = output.L_FA_Num+1;
                    output.L_FA_RT = [output.L_FA_RT, (FA_RT+250+double(input.stimOnTimeMs))];
                end
                
                
            end
       
        
    end
end
% calculate CR number, get rid of those followed by target

output.CR_Num = 0;
output.S_CR_Num = 0;
output.L_CR_Num =0;

for i = 1:trialnum
    cycle = [];
   
    CR_offs = [];
    cycle_i = [];
    cycle =tCycleNum(i);
    
    
    % separate out early and non early conditions
    
    
    if cycle>=3
        % if it is early trial,
        if Early(i)
            Triallength =Leverup(i)-Leverdown(i);
            FA_RT = (Triallength-250.*(cycle-1)-(double(input.stimOnTimeMs)*(cycle-1)));
            % if it is super early release
            if FA_RT<200 
                CR_offs =length(2:(cycle-3));
                cycle_i = 2:(cycle-3);
                
                    output.CR_Num= output.CR_Num+CR_offs;
                    output.S_CR_Num= output.S_CR_Num+sum(cycle_i<=4);
                    output.L_CR_Num= output.L_CR_Num+sum(cycle_i>=6 & cycle_i<=8);
                    
              
            else
                CR_offs =length(2:(cycle-2));
                cycle_i = 2:(cycle-2);
               
                    output.CR_Num= output.CR_Num+CR_offs;
                    
                    output.S_CR_Num= output.S_CR_Num+sum(cycle_i<=4);
                    output.L_CR_Num= output.L_CR_Num+sum(cycle_i>=6 & cycle_i<=8);
               
                
                
            end
        end
        % if it is success or miss trial,CR is all the previouse cycles
        if ~Early(i)
            CR_offs =length(2:(cycle-1));
            cycle_i = 2:(cycle-1);
            
                output.CR_Num= output.CR_Num+CR_offs;
                output.S_CR_Num= output.S_CR_Num+sum(cycle_i<=4);
                output.L_CR_Num= output.L_CR_Num+sum(cycle_i>=6 & cycle_i<=8);
           
        end
    end
    
    
    
end
% calculate FA rate


    [c d]= binofit(output.FA_Num,(output.FA_Num+output.CR_Num));
    output.all.FA_confi(1,1:2) = d;
    output.all.FA(1,1)=c;
    
    [c d]= binofit(output.S_FA_Num,(output.S_FA_Num+output.S_CR_Num));
    output.S_FA_confi(1,1:2) = d;
    output.S_FA(1,1)=c;
    
    [c d]= binofit(output.L_FA_Num,(output.L_FA_Num+output.L_CR_Num));
    output.L_FA_confi(1,1:2) = d;
    output.L_FA(1,1)=c;
    







end