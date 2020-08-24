function [output] = ISI_Samp(input,Off,Orien,Success,Miss,Early,tCycleNum,tOrientation,trialnum,Leverup,Leverdown)
idx_early = find(Early==1);
idx_HM = find(Early==0);
Triallength = nan(1,trialnum);
% get the trial length for all hit and mises trials
for i = 1:length(idx_HM)
    idx_temp = [];
    cycle = [];
    trialoffs= [];
    idx_temp = idx_HM(i);% Success and miss trials
    cycle = tCycleNum(idx_temp);
    if length(Off)>1
        trialoffs = double(input.tStimOffTimes{idx_temp});
        Triallength(idx_temp) = sum(trialoffs(1:cycle)) + double(input.stimOnTimeMs)*cycle;
    else
        Triallength(idx_temp) = cycle.*250 + double(input.stimOnTimeMs)*cycle;
    end
end


% consider trial length of fixed interval is shorter
% use different bins here
edges = [];


edges (1:length(950:500:2450)) = 950:500:2450;
edges(length(edges)+1) =8000;

cbin = mean([edges(1:end-1); edges(2:end)]);
[~,ind] = histc(Triallength,edges);

output.cbin.Hit = cbin;
output.outcome = {};
% caclulate hit rate bined by trial length
for i_bin = 1:length(cbin)
    
    for i_orien = 1: length(Orien)
        % find all the corrects on orientation changes
        
        
        output.HT_num(i_bin,i_orien) = sum(Success & (ind==i_bin)&(tOrientation==Orien(i_orien)));
        output.Miss_num(i_bin,i_orien) = sum (Miss &(ind==i_bin)&(tOrientation==Orien(i_orien)));
        output.Trial_num(i_bin,i_orien) = output.HT_num(i_bin,i_orien)+ output.Miss_num(i_bin,i_orien);
        output.HT_rate(i_bin,i_orien)= output.HT_num(i_bin,i_orien) ./ (output.HT_num(i_bin,i_orien)+output.Miss_num(i_bin,i_orien));
        temp = [];
        temp(1:output.HT_num(i_bin,i_orien)) =  ones(1,output.HT_num(i_bin,i_orien));
        temp((output.HT_num(i_bin,i_orien)+1):output.Trial_num(i_bin,i_orien)) = zeros(1,output.Miss_num(i_bin,i_orien));
        temp_rand  = temp(randperm(length(temp)));
        output.outcome{i_bin,i_orien} = temp_rand;% leave this for match across led conditions as well
        
        
        
    end
    % throw away the orientations that has less than 5 trials in total
    % get the fit for the threshold
    if output.Trial_num(i_bin,1)>=5 && output.Trial_num(i_bin,end)>=5
        trialVec = [];
        trialVec = output.HT_num(i_bin,:)+output.Miss_num(i_bin,:);
        output.fit{i_bin} = weibullFitLG(Orien, output.HT_rate(i_bin,:),1, 1, {'nTrials', trialVec});
        output.thresh(i_bin) =output.fit{i_bin}.thresh;
    end 
    if output.Trial_num(i_bin,1)<5 && output.Trial_num(i_bin,end)>=5
        trialVec = [];
        trialVec = output.HT_num(i_bin,:)+output.Miss_num(i_bin,:);
        output.fit{i_bin} = weibullFitLG(Orien(2:end), output.HT_rate(i_bin,2:end),1, 1, {'nTrials', trialVec(2:end)});
        output.thresh(i_bin) =output.fit{i_bin}.thresh;
    end
    if output.Trial_num(i_bin,1)<5 && output.Trial_num(i_bin,end)<5
        trialVec = [];
        trialVec = output.HT_num(i_bin,:)+output.Miss_num(i_bin,:);
        output.fit{i_bin} = weibullFitLG(Orien(2:end-1), output.HT_rate(i_bin,2:end-1),1, 1, {'nTrials', trialVec(2:end-1)});
        output.thresh(i_bin) =output.fit{i_bin}.thresh;
    end 
    if output.Trial_num(i_bin,1)>=5 && output.Trial_num(i_bin,end)<5
        trialVec = [];
        trialVec = output.HT_num(i_bin,:)+output.Miss_num(i_bin,:);
        output.fit{i_bin} = weibullFitLG(Orien(1:end-1), output.HT_rate(i_bin,1:end-1),1, 1, {'nTrials', trialVec(1:end-1)});
        output.thresh(i_bin) =output.fit{i_bin}.thresh;
    end 
    
end



% caclulate FA lengh on different offs


LengthFA = [];


for i = 1: length(idx_early)
    idx_temp = [];
    cycle = [];
    trialoffs= [];
    Triallength = [];
    FA_RT = [];
    idx_temp = idx_early(i);% early trials
    cycle = tCycleNum(idx_temp);
    Triallength =Leverup(idx_temp)-Leverdown(idx_temp);
    trialoffs = [];
    if length(Off)>1
        trialoffs = double(input.tStimOffTimes{idx_temp});
    else
        trialoffs =  repmat(250,1,cycle);
    end
    FA_RT = (Triallength-sum(trialoffs(1:(cycle-1)))-(double(input.stimOnTimeMs)*(cycle-1)));
    
    if cycle>=3
        if FA_RT>=200 && FA_RT<=550
            temp = [];
            
            temp = sum(trialoffs(1:(cycle-1)))+double(input.stimOnTimeMs)*(cycle-1);
            
            
            LengthFA = [LengthFA,temp];
            
        end
        
        
        if FA_RT<200 && trialoffs(cycle-1)==250
            
            LengthFA = [LengthFA,(sum(trialoffs(1:(cycle-2)))+double(input.stimOnTimeMs)*(cycle-2))];
            
        end
        
        
    end
end
% calculate CR trial length on different offs


LengthCR = [];

for i = 1:trialnum
    cycle = [];
    trialoffs= [];
    CR_offs = [];
    cycle =tCycleNum(i);
    if length(Off)>1
        trialoffs = double(input.tStimOffTimes{i});
    else
        trialoffs = repmat(250,1,cycle);
    end
    % separate out early and non early conditions
    
    if cycle>=3
        % if it is early trial,
        if Early(i)
            Triallength =Leverup(i)-Leverdown(i);
            FA_RT = (Triallength-sum(trialoffs(1:(cycle-1)))-(double(input.stimOnTimeMs)*(cycle-1)));
            % if it is super early release
            if FA_RT<200 && trialoffs(cycle-1)==250
                CR_offs =trialoffs(2:(cycle-3));
                cycle_i = 2:(cycle-3);
                temp=zeros(1,length(CR_offs));
                if ~isempty(CR_offs)
                    for j = 1:length(CR_offs)
                        temp(j) = sum(trialoffs(1:cycle_i(j)))+double(input.stimOnTimeMs)*cycle_i(j);
                    end
                    LengthCR = [LengthCR, temp];
                end
            else
                CR_offs =trialoffs(2:(cycle-2));
                cycle_i = 2:(cycle-2);
                temp=zeros(1,length(CR_offs));
                for j = 1:length(CR_offs)
                    temp(j) = sum(trialoffs(1:cycle_i(j)))+double(input.stimOnTimeMs)*cycle_i(j);
                end
                LengthCR = [LengthCR, temp];
                
                
                
            end
        end
        % if it is success or miss trial,CR is all the previouse cycles
        if ~Early(i)
            CR_offs =trialoffs(2:(cycle-1));
            cycle_i = 2:(cycle-1);
            temp=zeros(1,length(CR_offs));
            for j = 1:length(CR_offs)
                temp(j) = sum(trialoffs(1:cycle_i(j)))+double(input.stimOnTimeMs)*cycle_i(j);
            end
            LengthCR = [LengthCR, temp];
            
        end
    end
    
    
    
end

% cacllate the FA rate based on bines
edges = [];

edges (1:length(950:500:2450)) = 950:500:2450;
edges(length(edges)+1) =8000;

cbin = mean([edges(1:end-1); edges(2:end)]);
[CR_bin,ind_CR] = histc(LengthCR,edges);
[FA_bin,ind_FA] = histc(LengthFA,edges);
output.FA = FA_bin(1:end-1)./(FA_bin(1:end-1)+CR_bin(1:end-1));
output.cbin.FA = cbin;




end