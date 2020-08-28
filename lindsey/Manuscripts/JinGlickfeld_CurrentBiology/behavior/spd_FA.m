function [output] = spd_FA(Early,trialnum,HoldTimeMs,boundary)

idx_early = find(Early==1);
output.FA_Num = 0;
output.S_FA_Num = 0;
output.L_FA_Num = 0;
output.FA_RT = [];
for i = 1: length(idx_early)
    idx_temp = [];
 
    
    Triallength = [];
    FA_RT = [];
    idx_temp = idx_early(i);% early trials
    
    Triallength =HoldTimeMs(idx_temp);
    RT = mod(Triallength,500);
    
    if Triallength>=boundary  % only consider FA occur after a true target would occur
        
        % no RT for FA trial for now...        
        output.FA_Num = output.FA_Num+1;
        
        if  Triallength<=2700   % 500*4+700 cycle number 4 short trials
            output.S_FA_Num = output.S_FA_Num+1;            
        end
        if Triallength>=3200  % 500*6+200  long trials 
            output.L_FA_Num = output.L_FA_Num+1;            
        end
        
         if RT<200
           output.FA_RT = [output.FA_RT, RT+500];
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
    % separate out early and non early conditions
    RT = mod(Triallength,500); % 100+400 this is determined by the reaction time window for targets
    cycle = floor(Triallength./500);
    
    if Triallength>=boundary
        % if it is early trial,
        if Early(i)
           
            % if it is super early release
            if RT<200 
                CR_offs =length(1:(cycle-2));
                cycle_i = 1:(cycle-2);
                
                    output.CR_Num= output.CR_Num+CR_offs;
                    output.S_CR_Num= output.S_CR_Num+sum(cycle_i<=3);
                    output.L_CR_Num= output.L_CR_Num+sum(cycle_i>=5);
                    
              
            else
                CR_offs =length(1:(cycle-1));
                cycle_i = 1:(cycle-1);
               
                    output.CR_Num= output.CR_Num+CR_offs;
                    
                    output.S_CR_Num= output.S_CR_Num+sum(cycle_i<=3);
                    output.L_CR_Num= output.L_CR_Num+sum(cycle_i>=5);
               
                
                
            end
        end
        % if it is success or miss trial,CR is all the previouse cycles
        if ~Early(i)
            
            
            CR_offs =length(1:(cycle-2));
            cycle_i = 1:(cycle-2);
            
                output.CR_Num= output.CR_Num+CR_offs;
                output.S_CR_Num= output.S_CR_Num+sum(cycle_i<=3);
                output.L_CR_Num= output.L_CR_Num+sum(cycle_i>=5);
           
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