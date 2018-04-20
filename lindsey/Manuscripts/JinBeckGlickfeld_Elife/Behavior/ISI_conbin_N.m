function [output] = ISI_conbin_N(Early,input,Off,tCycleNum,trialnum,Leverup,Leverdown,Success,tOrientation,Orien,RT_HM,StimOff,StimOff_pre,Miss)

idx_early = find(Early==1);


% first off is off before the BS before the target, second off is the off
% before target.


output.HT_num = zeros(length(Off),length(Off),length(Orien));
output.Miss_num = zeros(length(Off),length(Off),length(Orien));
output.HT_rate = zeros(length(Off),length(Off),length(Orien));

for i_off = 1:length(Off) % the off before before
    for ii_off = 1:length(Off) % the off before the target
        for i_orien = 1: length(Orien)
            output.HT_num(i_off,ii_off,i_orien) = sum(Success & (StimOff==Off(ii_off))&(StimOff_pre==Off(i_off))&(tOrientation==Orien(i_orien)));
            output.Miss_num(i_off,ii_off,i_orien) = sum (Miss &(StimOff==Off(ii_off))&(StimOff_pre==Off(i_off))&(tOrientation==Orien(i_orien)));
            output.HT_rate(i_off,ii_off,i_orien) = output.HT_num(i_off,ii_off,i_orien) ./ (output.HT_num(i_off,ii_off,i_orien)+output.Miss_num(i_off,ii_off,i_orien));
            
            [a,b]= binofit(output.HT_num(i_off,ii_off,i_orien),(output.HT_num(i_off,ii_off,i_orien)+output.Miss_num(i_off,ii_off,i_orien)));
            output.confi(i_off,ii_off,i_orien,1:2) = b;
            output.c_hit(i_off,ii_off,i_orien) = a;
            
            output.RTonHit{i_off,1}{ii_off,1}{i_orien,1} = RT_HM(Success &(StimOff==Off(ii_off))&(StimOff_pre==Off(i_off))&(tOrientation==Orien(i_orien)));
            output.RTonHit_mean {i_off,1}{ii_off,1}(1,i_orien) = mean(output.RTonHit {i_off,1}{ii_off,1}{i_orien,1});
            output.RTonHit_mean {i_off,1}{ii_off,1}(2,i_orien) = std(output.RTonHit {i_off,1}{ii_off,1}{i_orien,1})./sqrt(output.HT_num(i_off,ii_off,i_orien));
            
        end
    end
end


% get the FA rate from different offs

for i_off = 1:length(Off)
    for ii_off = 1:length(Off)
        output.FA_RT{i_off,ii_off}=[];
    end
end
output.FA_Num = zeros(length(Off),length(Off));


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
    
    for  i_off = 1:length(Off) % off pre pre
        if cycle>=3 && cycle==length(trialoffs)  % only consider FA occur after cycle 2, because target never occur on cycle 1&2; through out early release for target
            %calculate collapsed FAs
            
            for ii_off = 1:length(Off) % off before the last BS
                % if release in the time window-extended to the 200 - 550
                if trialoffs(cycle-1)==Off(ii_off) && FA_RT>=200 && FA_RT<=550
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
                    if temp==Off(i_off)
                        output.FA_Num(i_off,ii_off) = output.FA_Num(i_off,ii_off)+1;
                        output.FA_RT{i_off,ii_off} = [output.FA_RT{i_off,ii_off}, FA_RT];
                        
                    end
                    
                    
                end
                % if release very fast and previous off time is 250ms off, then
                % this is count as FA for the previous cycle..
                
                if trialoffs(cycle-2)==Off(ii_off) && FA_RT<200 && trialoffs(cycle-1)==250 && cycle>3
                    temp = 0; % throw away trials that is randomized
                    if length(unique(trialoffs(1:(cycle-3)))) ==1
                        if unique(trialoffs(1:(cycle-3)))==250
                            temp=250;
                        end
                        if unique(trialoffs(1:(cycle-3)))==750
                            temp=750;
                        end
                        if unique(trialoffs(1:(cycle-3)))==500
                            temp=500;
                        end
                    end
                    if  temp==Off(i_off) % note that we lose some cycle 3 trials
                        
                        
                        output.FA_Num(i_off,ii_off) = output.FA_Num(i_off,ii_off)+1;
                        output.FA_RT{i_off,ii_off}  = [output.FA_RT{i_off,ii_off} , (FA_RT+250+double(input.stimOnTimeMs))];
                    end
                end
                
                
                
            end
        end
        
    end
end

% calculate CR number, get rid of those followed by target

output.CR_Num = zeros(length(Off),length(Off));



for i = 1:trialnum
    cycle = [];
    trialoffs= [];
    CR_offs = [];
    CR_pres=[];
    CR_presbin = [];
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
                if cycle>4
                    
                    
                    
                    CR_offs =trialoffs(2:(cycle-3));
                    % only modify the pres
                    
                    CR_pres = trialoffs(1:(cycle-4));
                    CR_presbin = [];
                    for i_cycle = 1:length(CR_pres)
                        temp = 0;
                        if length(unique(CR_pres(1:i_cycle))) ==1
                        if unique(CR_pres(1:i_cycle))==250
                            temp=250;
                        end
                        if unique(CR_pres(1:i_cycle))==500
                            temp=500;
                        end
                        if unique(CR_pres(1:i_cycle))==750
                            temp=750;
                        end
                        
                        end
                      
                        
                        CR_presbin(i_cycle) = temp;
                    end
                    
                    for i_off = 1:length(Off) % pre pre off
                        for ii_off =  1:length(Off) % pre off
                            output.CR_Num(i_off,ii_off)= output.CR_Num(i_off,ii_off)+sum(CR_presbin==Off(i_off)& CR_offs==Off(ii_off));
                        end
                        
                    end
                end
                
            else
                
                
                
                CR_offs =trialoffs(2:(cycle-2));
                
                
                CR_pres = trialoffs(1:(cycle-3));
                CR_presbin = [];
                if length(CR_pres)~=0
                    
                    
                    for i_cycle = 1:length(CR_pres)
                        
                         temp = 0;
                        if length(unique(CR_pres(1:i_cycle))) ==1
                        if unique(CR_pres(1:i_cycle))==250
                            temp=250;
                        end
                        if unique(CR_pres(1:i_cycle))==500
                            temp=500;
                        end
                        if unique(CR_pres(1:i_cycle))==750
                            temp=750;
                        end
                        
                        end
                      
                        
                       
                        CR_presbin(i_cycle) = temp;
                    end
                else
                    CR_presbin =CR_pres;
                end
                
                
                for i_off = 1:length(Off) % pre pre off
                    for ii_off =  1:length(Off) % pre off
                        output.CR_Num(i_off,ii_off)= output.CR_Num(i_off,ii_off)+sum(CR_presbin==Off(i_off)& CR_offs==Off(ii_off));
                        
                    end
                end
                
                
            end
        end
        % if it is success or miss trial,CR is all the previouse cycles
        if ~Early(i)
            CR_offs =trialoffs(2:(cycle-1));
            CR_pres = trialoffs(1:(cycle-2));
            CR_presbin = [];
            if length(CR_pres)~=0
                for i_cycle = 1:length(CR_pres)                    
                    temp = 0;
                    if length(unique(CR_pres(1:i_cycle))) ==1
                        if unique(CR_pres(1:i_cycle))==250
                            temp=250;
                        end
                        if unique(CR_pres(1:i_cycle))==500
                            temp=500;
                        end
                        if unique(CR_pres(1:i_cycle))==750
                            temp=750;
                        end
                        
                    end
                    
                    
                    CR_presbin(i_cycle) = temp;
                end
            else
                CR_presbin =CR_pres;
            end
            
            for i_off = 1:length(Off)
                for ii_off =  1:length(Off) % pre off
                    output.CR_Num(i_off,ii_off)= output.CR_Num(i_off,ii_off)+sum(CR_presbin==Off(i_off)& CR_offs==Off(ii_off));
                end
            end
        end
    end
    
    
    
end
% calculate FA rate

for i_off = 1:length(Off)
    for ii_off = 1:length(Off)
        [c,d]= binofit(output.FA_Num(i_off,ii_off),(output.FA_Num(i_off,ii_off)+output.CR_Num(i_off,ii_off)));
        output.FA_confi(i_off,ii_off,1:2) = d;
        output.FA(i_off,ii_off)=c;
        
        
    end
end

% get the FA rate for different conditions
temp =[];
temp1 =[];
temp = sum(output.FA_Num,1);
temp1 = sum(output.CR_Num,1);

for i_off = 1:length(Off)
    [a,b] = binofit(temp(i_off),temp(i_off)+temp1(i_off));
    output.pre.FA(i_off,1) = a;
    output.pre.FA_confi(i_off,1:2) = b;
    output.pre.FA_RT{i_off,1} = cat(2,output.FA_RT{:,i_off});
end
temp=[];
temp1=[];
temp = sum(output.FA_Num,2);
temp1 = sum(output.CR_Num,2);

for i_off = 1:length(Off)
    [a,b] = binofit(temp(i_off),temp(i_off)+temp1(i_off));
    output.pre2.FA(i_off,1) = a;
    output.pre2.FA_confi(i_off,1:2) = b;
    output.pre2.FA_RT{i_off,1} = cat(2,output.FA_RT{i_off,:});
end









end