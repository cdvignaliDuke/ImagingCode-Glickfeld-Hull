function [output] = ISI_RT(Early,input,Off,tCycleNum,Leverup,Leverdown,Success,Orien,tOrientation)

idx_early = find(Early==1);
idx_success = find(Success==1);

Orien_idx = zeros(1,length(tOrientation));
for i_orien = 1:length(Orien)
Orien_idx(tOrientation==(Orien(i_orien))) = i_orien;
output.T_orien {i_orien,1} = [];
end 

for i_off = 1:length(Off)
    output.T_RT{i_off,1}=[];
    output.B_RT{i_off,1} = [];
    output.pre_eB_RT{i_off,1} = [];
    output.pre_B_RT{i_off,1} = [];
   
end


% consider early relsease first
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
    if cycle>=3 && cycle==length(trialoffs)  % only consider FA occur after cycle 2, because target never occur on cycle 1&2;
        
        for i_off = 1:length(Off)
            if trialoffs(cycle-1)==Off(i_off) && FA_RT<=550
                output.B_RT{i_off,1} = [output.B_RT{i_off,1}, FA_RT];
            end
        end
        
%         if trialoffs(cycle-2)==Off(i_off) && FA_RT<200 && trialoffs(cycle-1)==250
%             output.B_RT{i_off,1} = [output.B_RT{i_off,1}, FA_RT+350];
%         end
    end
    % for early release of previous BS
    if cycle>=4 && cycle==length(trialoffs)
        
        for i_off = 1:length(Off)
            if trialoffs(cycle-2)==Off(i_off)
                output.pre_eB_RT{i_off,1} = [output.pre_eB_RT{i_off,1},(FA_RT+trialoffs(cycle-1)+double(input.stimOnTimeMs))];
                
            end
        end
        
    end
    % for the early release of the target response
    if cycle < length(trialoffs)
        new_RT = Triallength-sum(trialoffs(1:cycle))-(double(input.stimOnTimeMs)*cycle);
        
        output.T_orien{Orien_idx(idx_temp),1} = [output.T_orien{Orien_idx(idx_temp),1},new_RT];
        for i_off = 1:length(Off)
            if trialoffs(cycle)==Off(i_off)
                
                output.T_RT{i_off,1} = [output.T_RT{i_off,1},new_RT];
                
                
            end
        end
    end
    
    
end
% calculate target response seperate  for all the orienatation changes

for i = 1: length(idx_success)
    idx_temp = [];
    cycle = [];
    trialoffs= [];
    Triallength = [];
    RT = [];
    idx_temp = idx_success(i);% early trials
    cycle = tCycleNum(idx_temp);
    Triallength =Leverup(idx_temp)-Leverdown(idx_temp);
    trialoffs = double(input.tStimOffTimes{idx_temp});
    RT = (Triallength-sum(trialoffs(1:cycle))-(double(input.stimOnTimeMs)*cycle));
   
    
    output.T_orien{Orien_idx(idx_temp),1} = [output.T_orien{Orien_idx(idx_temp),1}, RT];
    
   
        
    
end
% calculate target response concatenate by orientations

for i = 1: length(idx_success)
    idx_temp = [];
    cycle = [];
    trialoffs= [];
    Triallength = [];
    RT = [];
    idx_temp = idx_success(i);% early trials
    cycle = tCycleNum(idx_temp);
    Triallength =Leverup(idx_temp)-Leverdown(idx_temp);
    trialoffs = double(input.tStimOffTimes{idx_temp});
    RT = (Triallength-sum(trialoffs(1:cycle))-(double(input.stimOnTimeMs)*cycle));
   
    
    for i_off = 1:length(Off)
        if trialoffs(cycle)==Off(i_off)
           output.T_RT{i_off,1}=[output.T_RT{i_off,1},RT];
        end
    end
    
    if cycle>=3  
        
        for i_off = 1:length(Off)
            
            if trialoffs(cycle-1)==Off(i_off)
                output.pre_B_RT{i_off,1} = [output.pre_B_RT{i_off,1}, (RT+trialoffs(cycle)+double(input.stimOnTimeMs))];
            end
        end
    end
   
        
    
end



for i_off = 1:length(Off)
     output.pre_B_RT{i_off,1} = [output.pre_B_RT{i_off,1},output.pre_eB_RT{i_off,1}];
end


end