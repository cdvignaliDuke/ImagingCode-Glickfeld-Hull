function [output] = ISI_Trl(StimOff,input,Off,Orien,Success,Miss,Early,tCycleNum,tOrientation,trialnum,Leverup,Leverdown)
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
    
    trialoffs = double(input.tStimOffTimes{idx_temp});
    Triallength(idx_temp) = sum(trialoffs(1:cycle)) + double(input.stimOnTimeMs)*cycle;
    
end
output.HM_length = Triallength; 


% get the mean trial length for each offs, get rid of nans 

for i_off = 1:length(Off) 
    
    output.length.HM{i_off,1} = Triallength(StimOff==Off(i_off));
    output.length.Hit(i_off) = nanmean(Triallength(StimOff==Off(i_off)));
    output.ori{i_off,1} = tOrientation(StimOff==Off(i_off));
    output.success{i_off,1} = Success(StimOff==Off(i_off));
    output.miss{i_off,1} = Miss(StimOff==Off(i_off));
end



% bin the trial length

edges = [];
edges(1) = 700;

% edges (2:(length(1200:300:5700)+1)) = 1200:300:5700;
edges (2:(length(1200:500:5700)+1)) = 1200:500:5700;
edges(length(edges)+1) =8000;

cbin = mean([edges(1:end-1); edges(2:end)]);
[bincounts,ind] = histc(Triallength,edges);

output.cbin.Hit = cbin;
% caclulate hit rate bined by trial length
for i_bin = 1:length(cbin)
    
    for i_orien = 1: length(Orien)
        % find all the corrects on orientation changes
        
        
        output.HT_num(i_bin,i_orien) = sum(Success & (ind==i_bin)&(tOrientation==Orien(i_orien)));
        output.Miss_num(i_bin,i_orien) = sum (Miss &(ind==i_bin)&(tOrientation==Orien(i_orien)));
        output.HT_rate(i_bin,i_orien)= output.HT_num(i_bin,i_orien) ./ (output.HT_num(i_bin,i_orien)+output.Miss_num(i_bin,i_orien));
        
        
        
    end
    
    % get the fit for the threshold
    trialVec = [];
    trialVec = output.HT_num(i_bin,:)+output.Miss_num(i_bin,:);
    try
    output.fit{i_bin} = weibullFitLG(Orien, output.HT_rate(i_bin,:),1, 1, {'nTrials', trialVec});
    output.thresh(i_bin) =output.fit{i_bin}.thresh;
    catch
    output.thresh(i_bin) = NaN;
    end 
    
end




% caclulate FA lengh on different offs


for i_off = 1:length(Off)
    output.FA_length{i_off,1}=[];
    
end

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
    trialoffs = double(input.tStimOffTimes{idx_temp});
    FA_RT = (Triallength-sum(trialoffs(1:(cycle-1)))-(double(input.stimOnTimeMs)*(cycle-1)));
    if cycle>=3 && cycle==length(trialoffs)  % only consider FA occur after cycle 2, because target never occur on cycle 1&2; through out early release for target
        %calculate collapsed FAs
        
        for i_off = 1:length(Off)
            % if release in the RT window 200-550 ms
            if trialoffs(cycle-1)==Off(i_off) && FA_RT>=200 && FA_RT<=550
                LengthFA = [LengthFA,(sum(trialoffs(1:(cycle-1)))+double(input.stimOnTimeMs)*(cycle-1))];
                
                output.FA_length{i_off,1} = [output.FA_length{i_off,1},(sum(trialoffs(1:(cycle-1)))+double(input.stimOnTimeMs)*(cycle-1))];
                
            end
            % if release very fast and previous off time is 250ms off, then
            % this is count as FA for the previous cycle..
            
            if trialoffs(cycle-2)==Off(i_off) && FA_RT<200 && trialoffs(cycle-1)==250 && cycle>3
                
                LengthFA = [LengthFA,(sum(trialoffs(1:(cycle-2)))+double(input.stimOnTimeMs)*(cycle-2))];
                
                output.FA_length{i_off,1} = [output.FA_length{i_off,1}, (sum(trialoffs(1:(cycle-2)))+double(input.stimOnTimeMs)*(cycle-2))];
                
            end
        end
        
    end
end
% calculate CR trial length on different offs

for i_off = 1:length(Off)
    output.CR_length{i_off,1}=[];
    
end


LengthCR = [];

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
                cycle_i = 2:(cycle-3);
                temp=zeros(1,length(CR_offs));
                if ~isempty(CR_offs)
                for j = 1:length(CR_offs)
                    temp(j) = sum(trialoffs(1:cycle_i(j)))+double(input.stimOnTimeMs)*cycle_i(j);
                end
                LengthCR = [LengthCR, temp];
                for i_off = 1:length(Off)
                    
                    output.CR_length{i_off,1} = [output.CR_length{i_off,1},temp(CR_offs==Off(i_off))];
                    
                    
                end
                end
            else
                CR_offs =trialoffs(2:(cycle-2));
                cycle_i = 2:(cycle-2);
                temp=zeros(1,length(CR_offs));
                for j = 1:length(CR_offs)
                    temp(j) = sum(trialoffs(1:cycle_i(j)))+double(input.stimOnTimeMs)*cycle_i(j);
                end
                LengthCR = [LengthCR, temp];
                for i_off = 1:length(Off)
                    
                    output.CR_length{i_off,1} = [output.CR_length{i_off,1},temp(CR_offs==Off(i_off))];
                    
                    
                end
                
                
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
            for i_off = 1:length(Off)
                
                output.CR_length{i_off,1} = [output.CR_length{i_off,1},temp(CR_offs==Off(i_off))];
                
                
            end
        end
    end
    
    
    
end
% calculate the mean hold time for FA measurement for each offs
for i_off = 1:length(Off)
    output.length.FA(i_off) = mean([output.CR_length{i_off,1}, output.FA_length{i_off,1}]);
    
   
end 

% cacllate the FA rate based on bines


% bin the trial length
edges = [];
edges(1) = 700;
% edges (2:(length(1200:300:4500)+1)) = 1200:300:4500;
edges (2:(length(1200:500:4200)+1)) = 1200:500:4200;
edges(length(edges)+1) = 7000;
cbin = mean([edges(1:end-1); edges(2:end)]);
[CR_bin,ind_CR] = histc(LengthCR,edges);
[FA_bin,ind_FA] = histc(LengthFA,edges);
output.FA = FA_bin(1:end-1)./(FA_bin(1:end-1)+CR_bin(1:end-1));
output.cbin.FA = cbin;

%% calculate the hit and FA rate for three intervals in the same range
% first hit rate, range should be bigger than 1200 and less than 6200 %

edges = [];
edges (1:(length(1200:500:6200))) = 1200:500:6200;
cbin = [];
cbin = mean([edges(1:end-1); edges(2:end)]);

num_off = [];
% from each bin, get equal amount of trials within different off conditions
for i_off = 1:length(Off)
    [num_off(i_off,:),idxOff{i_off,1}] = histc(output.length.HM{i_off,1},edges);
    
end 
resample_num = min(num_off,[],1);

success_temp = {};
miss_temp = {};
orientation_temp = {};
for i_off = 1:length(Off)
    for i_bin = 1:length(cbin)
        idx_temp = [];
        idx_temp = find(idxOff{i_off,1}==i_bin);
        
        
        success_temp{i_off,i_bin} = output.success{i_off,1}(idx_temp(1:resample_num(i_bin))) ;
        
        miss_temp{i_off,i_bin} =  output.miss{i_off,1}(idx_temp(1:resample_num(i_bin))) ;
       
        orientation_temp{i_off,i_bin} = output.ori{i_off,1}(idx_temp(1:resample_num(i_bin))) ;
        output.new.HM_length{i_off,i_bin} = output.length.HM{i_off,1}(idx_temp(1:resample_num(i_bin)));
       
    end
end
% cat the trial of all the bins
for i_off = 1:length(Off)
    N_success = [];
    N_success = cat(2,success_temp{i_off,:});
    N_miss = [];
    N_miss =  cat(2,miss_temp{i_off,:});
    N_ori = [];
    N_ori = cat(2,orientation_temp{i_off,:});
    N_length = [];
    N_length = cat(2,output.new.HM_length{i_off,:});
    output.new.N_length{i_off,1} =  N_length ;
    
    output.new.HM_length_mean(i_off) = mean(N_length) ;
    for i_orien = 1:length(Orien)
        output.new.HT_num(i_off,i_orien) = sum(N_success &(N_ori==Orien(i_orien)));
        output.new.Miss_num(i_off,i_orien) = sum(N_miss &(N_ori==Orien(i_orien)));
        output.new.HT_rate(i_off,i_orien)=  output.new.HT_num(i_off,i_orien) ./ (output.new.HT_num(i_off,i_orien)+output.new.Miss_num(i_off,i_orien));
        
    end 
    trialVec = [];
    trialVec = output.new.HT_num(i_off,:)+output.new.Miss_num(i_off,:);
    output.new.fit{i_off} = weibullFitLG(Orien, output.new.HT_rate(i_off,:),1, 1, {'nTrials', trialVec});
    output.new.thresh(i_off) =output.new.fit{i_off}.thresh;
   
end 

%% get the matched trial length for FAs
output.CR_length
output.FA_length
edges = [];
edges (1:(length(1200:500:4200))) = 1200:500:4200;
cbin = [];
cbin = mean([edges(1:end-1); edges(2:end)]);


% from each bin, get equal amount of trials within different off conditions
for i_off = 1:length(Off)
    % change the CR and FA's into zeros and ones 
    
    [CR_num(i_off,:),CR_idx{i_off,1}] = histc(output.CR_length{i_off,1},edges);
    [FA_num(i_off,:),FA_idx{i_off,1}] = histc(output.FA_length{i_off,1},edges);
    
end 

FACR_num = CR_num+FA_num;

resample_num = min(FACR_num,[],1);

for i_off = 1:length(Off)
    for i_bin = 1:length(cbin)
        idx_temp = [];
        idx_temp = find(CR_idx{i_off,1}==i_bin);
        CRslength{i_off,i_bin} = output.CR_length{i_off,1}(idx_temp);
        CRs{i_off,i_bin} =zeros(1,length(idx_temp));
        idx_temp = [];
        idx_temp = find(FA_idx{i_off,1}==i_bin);
        FAslength{i_off,i_bin} = output.FA_length{i_off,1}(idx_temp);       
        FAs{i_off,i_bin} = ones(1,length(idx_temp));
        % cat the data
        FACRslength{i_off,i_bin} = [CRslength{i_off,i_bin} FAslength{i_off,i_bin}];
        FACRs{i_off,i_bin} = [CRs{i_off,i_bin},FAs{i_off,i_bin} ];
   
    end
end

for i_off = 1:length(Off)
    for i_bin = 1:length(cbin)
        temp_idx = [];
        [N_FACRs{i_off,i_bin},temp_idx] = datasample(FACRs{i_off,i_bin},resample_num(i_bin),'Replace',false);
        N_FACRslength{i_off,i_bin} = FACRslength{i_off,i_bin} (temp_idx);
    end 
    NC_FRCRs{i_off,1} = cat(2,N_FACRs{i_off,:});
    NC_FRCRslength{i_off,1} = cat(2, N_FACRslength{i_off,:});
end 

% get the FA rate
for i_off = 1:length(Off)
output.new.CR_num(i_off) =  sum(~NC_FRCRs{i_off,1});
output.new.FA_num(i_off) =  sum(NC_FRCRs{i_off,1});
output.new.FA(i_off) = output.new.FA_num(i_off)/(output.new.CR_num(i_off) + output.new.FA_num(i_off));
output.new.FACR_lengthmean(i_off) = mean( NC_FRCRslength{i_off,1});
end

end