function [output] = ISI_FA_bin(Early,input,Off,tCycleNum,trialnum,Leverup,Leverdown)

idx_early = find(Early==1);
output.FA_Num = zeros(length(Off),1);

for i_off = 1:length(Off)
    output.FA_RT{i_off,1}=[];
   
    output.cycle{i_off,1} = [];
end


for i = 1: length(idx_early)
    idx_temp = [];
    cycle = [];
    trialoffs= [];
    Triallength = [];
    FA_RT = [];
    idx_temp = idx_early(i);% early trials
    cycle = tCycleNum(idx_temp);
    Triallength =Leverup(idx_temp)-Leverdown(idx_temp);
    trialoffs = double(input.tStimOffTimes{idx_temp});
    FA_RT = (Triallength-sum(trialoffs(1:(cycle-1)))-(double(input.stimOnTimeMs)*(cycle-1)));
    if cycle>=3 && cycle==length(trialoffs)&& FA_RT>=200 && FA_RT<=550  % only consider FA occur after cycle 2, because target never occur on cycle 1&2; through out early release for target
        %calculate collapsed FAs
        temp = 0; % throw away trials that is randomized
        if length(unique(trialoffs(1:(cycle-1)))) ==1
        if unique(trialoffs(1:(cycle-1)))==250
            temp=250;
        end
        if unique(trialoffs(1:(cycle-1)))==750
            temp=750;
        end
        if unique(trialoffs(1:(cycle-1)))==500
            temp=500;
        end  
        end
            
        for i_off = 1:length(Off)
            % if release in the RT window 200-550 ms           
            if temp==Off(i_off) 
                           
              
                output.FA_Num(i_off,1) = output.FA_Num(i_off,1)+1;
                output.FA_RT{i_off,1} = [output.FA_RT{i_off,1}, FA_RT];
                output.cycle{i_off,1} = [output.cycle{i_off,1}, cycle];
                
               
                
            end
        end 
    end 
    
    if cycle>=3 && cycle==length(trialoffs)&& FA_RT<200 && trialoffs(cycle-1)==250
        % if release very fast and previous off time is 250ms off, then
        % this is count as FA for the previous cycle..
        temp = 0; % throw away trials that is randomized
        if length(unique(trialoffs(1:(cycle-2)))) ==1
        if unique(trialoffs(1:(cycle-2)))==250
            temp=250;
        end
        if unique(trialoffs(1:(cycle-2)))==750
            temp=750;
        end
        if unique(trialoffs(1:(cycle-2)))==500
            temp=500;
        end  
        end
            
        
        for i_off = 1:length(Off)
            if temp==Off(i_off)
                output.FA_Num(i_off,1) = output.FA_Num(i_off,1)+1;
                output.FA_RT{i_off,1} = [output.FA_RT{i_off,1}, (FA_RT+250+double(input.stimOnTimeMs))];
                output.cycle{i_off,1} = [output.cycle{i_off,1}, cycle-1];
                
                
                
                
            end
        end
        
    end

end
% calculate CR number, get rid of those followed by target

output.CR_Num = zeros(length(Off),1);


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
            FA_RT = (Triallength-sum(trialoffs(1:(cycle-1)))-(double(input.stimOnTimeMs)*(cycle-1)));
            % if it is super early release
            if FA_RT<200 && trialoffs(cycle-1)==250
                CR_offs =trialoffs(2:(cycle-3));
                
                
                % need to the mean of previous CR_offs..
                
                for i_CR = 1:length(CR_offs)
                     
                    if length(unique(trialoffs(1:i_CR+1)))==1
                    
                    for i_off = 1:length(Off)
                    output.CR_Num(i_off,1)= output.CR_Num(i_off,1)+sum(unique(trialoffs(1:i_CR+1))==Off(i_off));
                    
                    
                    end
                    end 
                end 
            else
                CR_offs =trialoffs(2:(cycle-2));
               
                for i_CR = 1:length(CR_offs)
                if length(unique(trialoffs(1:i_CR+1)))==1
                for i_off = 1:length(Off)
                    output.CR_Num(i_off,1)= output.CR_Num(i_off,1)+sum(unique(trialoffs(1:i_CR+1))==Off(i_off));
                    
                    
                end
                end 
                end 
                
            end
        end
        % if it is success or miss trial,CR is all the previouse cycles
        if ~Early(i)
            CR_offs =trialoffs(2:(cycle-1));
            for i_CR = 1:length(CR_offs)
            if length(unique(trialoffs(1:i_CR+1)))==1
            
            for i_off = 1:length(Off)
                output.CR_Num(i_off,1)= output.CR_Num(i_off,1)+sum(unique(trialoffs(1:i_CR+1))==Off(i_off));
                
            end
            end
            end
        end
    end
    
    
    
end
% calculate FA rate

for i_off = 1:length(Off)
    [c d]= binofit(output.FA_Num(i_off,1),(output.FA_Num(i_off,1)+output.CR_Num(i_off,1)));
    output.FA_confi(i_off,1:2) = d;
    output.FA(i_off,1)=c;
    
  
    

    
end
% caclulate FA collapsed all conditions 

[c d]= binofit(output.FA_Num(1,1)+output.FA_Num(2,1)+output.FA_Num(3,1),(output.FA_Num(1,1)+output.CR_Num(1,1)+output.FA_Num(2,1)+output.CR_Num(2,1)+output.FA_Num(3,1)+output.CR_Num(3,1)));
output.all.FA_confi(1:2) = d;
output.all.FA=c;


% calculate d prime and criterion for different orientations
%Adjust only the extreme values by replacing rates of 0 with 0.5/n0.5/n and rates
%of 1 with (n?0.5)/n(n?0.5)/n where nn is the number of signal or noise trials (Macmillan & Kaplan, 1985)
% for i_off= 1:length(Off)
%     for i_orien = 1:length(Orien)
%         if output.c_hit{i_off,1}(i_orien,1)<1
%             [output.dprime{i_off,1}(i_orien,1),output.criterion{i_off,1}(i_orien,1)] = dprime_simple(output.c_hit{i_off,1}(i_orien,1),output.FA(i_off,1));
%         else
%             signaltrialN = output.HT_num{i_off,1}(i_orien,1) + output.Miss_num{i_off,1}(i_orien,1);
%             [output.dprime{i_off,1}(i_orien,1),output.criterion{i_off,1}(i_orien,1)] = dprime_simple((0.5-signaltrialN)./signaltrialN,output.FA(i_off,1));
%         end
%     end
% end






end