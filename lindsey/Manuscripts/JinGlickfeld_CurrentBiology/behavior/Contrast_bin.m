function [output] = Contrast_bin(input,Contrast,Success,Miss,Early,tContrast,trialnum,HoldTimeMs,ReactTime,boundary)

idx_early = find(Early==1);
idx_HM = find(Early==0);
Triallength = nan(1,trialnum);
% get the trial length for all hit and mises trials
for i = 1:length(idx_HM)
    idx_temp = [];    
    idx_temp = idx_HM(i);% Success and miss trials        
    Triallength(idx_temp) = HoldTimeMs(idx_temp)-ReactTime(idx_temp);
   
end


% consider trial length of fixed interval is shorter
% use different bins here
edges = [];

edges = [700 1750 2450 3150 5000]; 
% modify bins later if necessary 
% edges (1:length(950:500:2450)) = 950:500:2450;
% edges(length(edges)+1) =8000;
% edges (1:length(950:500:3450)) = 950:500:3450;
% edges(length(edges)+1) =5000;

cbin = mean([edges(1:end-1); edges(2:end)]);
[~,ind] = histc(Triallength,edges);

output.cbin.Hit = cbin;
output.outcome = {};
% caclulate hit rate bined by trial length
for i_bin = 1:length(cbin)
    
    for i_contrast = 1: length(Contrast)
        % find all the corrects on orientation changes
        
        
        output.HT_num(i_bin,i_contrast) = sum(Success & (ind==i_bin)&(tContrast==Contrast(i_contrast)));
        output.Miss_num(i_bin,i_contrast) = sum (Miss &(ind==i_bin)&(tContrast==Contrast(i_contrast)));
        output.Trial_num(i_bin,i_contrast) = output.HT_num(i_bin,i_contrast)+ output.Miss_num(i_bin,i_contrast);
        output.HT_rate(i_bin,i_contrast)= output.HT_num(i_bin,i_contrast) ./ (output.HT_num(i_bin,i_contrast)+output.Miss_num(i_bin,i_contrast));
        temp = [];
        temp(1:output.HT_num(i_bin,i_contrast)) =  ones(1,output.HT_num(i_bin,i_contrast));
        temp((output.HT_num(i_bin,i_contrast)+1):output.Trial_num(i_bin,i_contrast)) = zeros(1,output.Miss_num(i_bin,i_contrast));
        temp_rand  = temp(randperm(length(temp)));
        output.outcome{i_bin,i_contrast} = temp_rand;% leave this for match across led conditions as well
        
        
        
    end
    % throw away the orientations that has less than 5 trials in total
    % get the fit for the threshold
    if output.Trial_num(i_bin,1)>=5 && output.Trial_num(i_bin,end)>=5
        trialVec = [];
        trialVec = output.HT_num(i_bin,:)+output.Miss_num(i_bin,:);
        output.fit{i_bin} = weibullFitLG(Contrast, output.HT_rate(i_bin,:),1, 1, {'nTrials', trialVec});
        output.thresh(i_bin) =output.fit{i_bin}.thresh;
    end 
    if output.Trial_num(i_bin,1)<5 && output.Trial_num(i_bin,end)>=5
        trialVec = [];
        trialVec = output.HT_num(i_bin,:)+output.Miss_num(i_bin,:);
        output.fit{i_bin} = weibullFitLG(Contrast(2:end), output.HT_rate(i_bin,2:end),1, 1, {'nTrials', trialVec(2:end)});
        output.thresh(i_bin) =output.fit{i_bin}.thresh;
    end
    if output.Trial_num(i_bin,1)<5 && output.Trial_num(i_bin,end)<5
        trialVec = [];
        trialVec = output.HT_num(i_bin,:)+output.Miss_num(i_bin,:);
        output.fit{i_bin} = weibullFitLG(Contrast(2:end-1), output.HT_rate(i_bin,2:end-1),1, 1, {'nTrials', trialVec(2:end-1)});
        output.thresh(i_bin) =output.fit{i_bin}.thresh;
    end 
    if output.Trial_num(i_bin,1)>=5 && output.Trial_num(i_bin,end)<5
        trialVec = [];
        trialVec = output.HT_num(i_bin,:)+output.Miss_num(i_bin,:);
        output.fit{i_bin} = weibullFitLG(Contrast(1:end-1), output.HT_rate(i_bin,1:end-1),1, 1, {'nTrials', trialVec(1:end-1)});
        output.thresh(i_bin) =output.fit{i_bin}.thresh;
    end 
    
end



% caclulate FA lengh 


LengthFA = [];

for i = 1: length(idx_early)
    idx_temp = [];
    Triallength = [];
    idx_temp = idx_early(i);% early trials
    Triallength = HoldTimeMs(idx_temp);
    if Triallength >= boundary % release after the fixed re
        LengthFA = [LengthFA,Triallength];
    end
    
end

% calculate CR (Miss+Hit) length 


LengthCR = [];

for i = 1:trialnum
    Triallength = [];
    Triallength = HoldTimeMs(i);
    
    if Early(i) && (Triallength - 550) >= boundary % 550 is the same as RT window 
     
    
    
    LengthCR = [LengthCR, (Triallength - 550)];
    end 
    if ~Early(i)
    LengthCR = [LengthCR, Triallength];    
        
        
    end 
end


% cacllate the FA rate based on bines
edges = [];
% edges (1:length(950:500:3450)) = 950:500:3450;
% edges(length(edges)+1) =5000;
edges = [];

edges = [700 1750 2450 3150 5000]; 

% edges (1:length(950:500:2450)) = 950:500:2450;
% edges(length(edges)+1) =8000;

cbin = mean([edges(1:end-1); edges(2:end)]);
[CR_bin,ind_CR] = histc(LengthCR,edges);
for i_bin = 1:length(CR_bin)
    CR_bin_new(i_bin) = sum(CR_bin(i_bin:end));
end 

[FA_bin,ind_FA] = histc(LengthFA,edges);
output.FA = FA_bin(1:end-1)./(FA_bin(1:end-1)+CR_bin_new(1:end-1));
output.cbin.FA = cbin;




end