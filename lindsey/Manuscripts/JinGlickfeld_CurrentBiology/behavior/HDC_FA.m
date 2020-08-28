function [output] = HDC_FA(Early,trialnum,HoldTimeMs,boundary)

idx_early = find(Early==1);
output.FA_Num = 0;
output.S_FA_Num = 0;
output.L_FA_Num = 0;
output.FA_RT=[];

for i = 1: length(idx_early)
    idx_temp = [];
 
    
    Triallength = [];
    RT = [];
    idx_temp = idx_early(i);% early trials
    
    Triallength =HoldTimeMs(idx_temp);
    RT = mod(Triallength,350);
    
    if Triallength>=boundary  % only consider FA occur after cycle 2, because target never occur on cycle 1&2; through out early release for target
        %calculate collapsed FAs
        
        % no RT for FA trial for now...
        output.FA_Num = output.FA_Num+1;
        
        if  Triallength<=2300   % 350*5+550cycle number 3,4,5 as short cycles
            output.S_FA_Num = output.S_FA_Num+1;
        end
        if Triallength>=2650 && Triallength<=3700  %7*350+200 350*9+550cycle number 7,8,9 as long cycles
            output.L_FA_Num = output.L_FA_Num+1;
        end
        
        if RT<200
           output.FA_RT = [output.FA_RT, RT+350];
        else
            output.FA_RT = [output.FA_RT, RT];
        end
        
    end
    
end

% calculate CR number, get rid of those followed by target

output.CR_Num = 0;
output.S_CR_Num = 0;
output.L_CR_Num =0;

for i = 1:trialnum
       
    Triallength =HoldTimeMs(i);
    RT = [];
    % separate out early and non early conditions
    RT = mod(Triallength,350);
    cycle = floor(Triallength./350);
    
    if Triallength>=boundary
        % if it is early trial,
        if Early(i)
           
            % if it is super early release
            if RT<200 
                CR_offs =length(1:(cycle-2));
                cycle_i = 1:(cycle-2);
                
                    output.CR_Num= output.CR_Num+CR_offs;
                    output.S_CR_Num= output.S_CR_Num+sum(cycle_i<=4);
                    output.L_CR_Num= output.L_CR_Num+sum(cycle_i>=6 & cycle_i<=8);
                    
              
            else
                CR_offs =length(1:(cycle-1));
                cycle_i = 1:(cycle-1);
               
                    output.CR_Num= output.CR_Num+CR_offs;
                    
                    output.S_CR_Num= output.S_CR_Num+sum(cycle_i<=4);
                    output.L_CR_Num= output.L_CR_Num+sum(cycle_i>=6 & cycle_i<=8);
               
                
                
            end
        end
        % if it is success or miss trial,CR is all the previouse cycles
        if ~Early(i)
            
            
            CR_offs =length(1:(cycle-2));
            cycle_i = 1:(cycle-2);
            
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