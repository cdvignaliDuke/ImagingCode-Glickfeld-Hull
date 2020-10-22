function [output] = ISI_Re(Early,input,Off,tCycleNum,trialnum,Leverup,Leverdown)

idx_early = find(Early==1);
output.Re_Num = zeros(length(Off),1);
output.S_Re_Num = zeros(length(Off),1);
output.L_Re_Num = zeros(length(Off),1);

for i_off = 1:length(Off)
    output.Re_RT{i_off,1}=[];
    output.S_Re_RT{i_off,1} = [];
    output.L_Re_RT{i_off,1} = [];
end


for i = 1: length(idx_early)
    idx_temp = [];
    cycle = [];
    trialoffs= [];
    Triallength = [];
    Re_RT = [];
    idx_temp = idx_early(i);% early trials
    cycle = tCycleNum(idx_temp);
    Triallength =Leverup(idx_temp)-Leverdown(idx_temp);
    trialoffs = double(input.tStimOffTimes{idx_temp});
    Re_RT = (Triallength-sum(trialoffs(1:(cycle-1)))-(double(input.stimOnTimeMs)*(cycle-1)));
    if cycle>=3 && cycle==length(trialoffs)  % only consider Re occur after cycle 2, because target never occur on cycle 1&2; through out early release for target
        %calculate collapsed Res
        
        for i_off = 1:length(Off)
            % if release during the 100 ms visual presentation 
            % need to exclude the true FAs 
            if trialoffs(cycle-1)==Off(i_off) && Re_RT>=0 && Re_RT<=100
              
                output.Re_Num(i_off,1) = output.Re_Num(i_off,1)+1;
                output.Re_RT{i_off,1} = [output.Re_RT{i_off,1}, Re_RT];
                
                if cycle<=5   % cycle number 3,4,5 as short cycles
                    output.S_Re_Num(i_off,1) = output.S_Re_Num(i_off,1)+1;
                    output.S_Re_RT{i_off,1} = [output.S_Re_RT{i_off,1}, Re_RT];
                end
                if cycle>=7 && cycle<=9  % cycle number 7,8,9 as long cycles
                    output.L_Re_Num(i_off,1) = output.L_Re_Num(i_off,1)+1;
                    output.L_Re_RT{i_off,1} = [output.L_Re_RT{i_off,1}, Re_RT];
                end
                
            end
            
%             % if release very fast and previous off time is 250ms off, then
%             % exclude this for 250ms off 
%             
%             if trialoffs(cycle-1)==Off(i_off) && Re_RT<50 && trialoffs(cycle-1)~=250
%                  output.Re_Num(i_off,1) = output.Re_Num(i_off,1)+1;
%                 output.Re_RT{i_off,1} = [output.Re_RT{i_off,1}, Re_RT];
%                 
%                 if cycle<=5   % cycle number 3,4,5 as short cycles
%                     output.S_Re_Num(i_off,1) = output.S_Re_Num(i_off,1)+1;
%                     output.S_Re_RT{i_off,1} = [output.S_Re_RT{i_off,1}, Re_RT];
%                 end
%                 if cycle>=7 && cycle<=9  % cycle number 7,8,9 as long cycles
%                     output.L_Re_Num(i_off,1) = output.L_Re_Num(i_off,1)+1;
%                     output.L_Re_RT{i_off,1} = [output.L_Re_RT{i_off,1}, Re_RT];
%                 end
%                 
%                 
%                 
%             end
            
            
            
          
        end
        
    end
end
% calculate CR number, get rid of those followed by target

output.CR_Num = zeros(length(Off),1);
output.S_CR_Num = zeros(length(Off),1);
output.L_CR_Num = zeros(length(Off),1);

for i = 1:trialnum
    cycle = [];
    trialoffs= [];
    CR_offs = [];
    cycle =tCycleNum(i);
    trialoffs = double(input.tStimOffTimes{i});
    
    % separate out early and non early conditions
    
    
    if cycle>=3
        % if it is early trial,
        if Early(i)
            Triallength =Leverup(i)-Leverdown(i);
            Re_RT = (Triallength-sum(trialoffs(1:(cycle-1)))-(double(input.stimOnTimeMs)*(cycle-1)));
            % if response is not during the visual presentations.
           
            if Re_RT>100
             CR_offs =trialoffs(2:(cycle-1));
            cycle_i = 2:(cycle-1);
            for i_off = 1:length(Off)
                output.CR_Num(i_off,1)= output.CR_Num(i_off,1)+sum(CR_offs==Off(i_off));
                output.S_CR_Num(i_off,1)= output.S_CR_Num(i_off,1)+sum(CR_offs(cycle_i<=4)==Off(i_off));
                output.L_CR_Num(i_off,1)= output.L_CR_Num(i_off,1)+sum(CR_offs(cycle_i>=6 & cycle_i<=8)==Off(i_off));
            end   
                
                
                

            else
                
                CR_offs =trialoffs(2:(cycle-2));
                cycle_i = 2:(cycle-2);
                for i_off = 1:length(Off)
                    output.CR_Num(i_off,1)= output.CR_Num(i_off,1)+sum(CR_offs==Off(i_off));
                    
                    output.S_CR_Num(i_off,1)= output.S_CR_Num(i_off,1)+sum(CR_offs(cycle_i<=4)==Off(i_off));
                    output.L_CR_Num(i_off,1)= output.L_CR_Num(i_off,1)+sum(CR_offs(cycle_i>=6 & cycle_i<=8)==Off(i_off));
                end
                
                
            end
        end
        % if it is success or miss trial,CR is all the previouse cycles
        if ~Early(i)
            CR_offs =trialoffs(2:(cycle-1));
            cycle_i = 2:(cycle-1);
            for i_off = 1:length(Off)
                output.CR_Num(i_off,1)= output.CR_Num(i_off,1)+sum(CR_offs==Off(i_off));
                output.S_CR_Num(i_off,1)= output.S_CR_Num(i_off,1)+sum(CR_offs(cycle_i<=4)==Off(i_off));
                output.L_CR_Num(i_off,1)= output.L_CR_Num(i_off,1)+sum(CR_offs(cycle_i>=6 & cycle_i<=8)==Off(i_off));
            end
        end
    end
    
    
    
end
% calculate Re rate

for i_off = 1:length(Off)
    [c d]= binofit(output.Re_Num(i_off,1),(output.Re_Num(i_off,1)+output.CR_Num(i_off,1)));
    output.Re_confi(i_off,1:2) = d;
    output.Re(i_off,1)=c;
    
    [c d]= binofit(output.S_Re_Num(i_off,1),(output.S_Re_Num(i_off,1)+output.S_CR_Num(i_off,1)));
    output.S_Re_confi(i_off,1:2) = d;
    output.S_Re(i_off,1)=c;
    
    [c d]= binofit(output.L_Re_Num(i_off,1),(output.L_Re_Num(i_off,1)+output.L_CR_Num(i_off,1)));
    output.L_Re_confi(i_off,1:2) = d;
    output.L_Re(i_off,1)=c;
    
end








end