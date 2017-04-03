function [trial_resp, trial_target, trial_length] = trialLengthMsMat(trL,tr,minTrLMs,trType,targets)

% find the trials that meet the cutoff (short or long only, or all trials)
if strcmp(trType,'short')
    ind = cellfun(@(x) x < minTrLMs, trL, 'unif',0);
    temp_tr = cellfun(@(x,y) x(:,:,y), tr, ind, 'unif',0);
elseif strcmp(trType,'long')
    ind = cellfun(@(x) x > minTrLMs, trL, 'unif',0);
    temp_tr = cellfun(@(x,y) x(:,:,y), tr, ind, 'unif',0);
elseif strcmp(trType,'all')
    ind = cellfun(@(x) x > 0, trL, 'unif',0);
    temp_tr = cellfun(@(x,y) x(:,:,y), tr, ind, 'unif',0);    
end

% fill in the new trial length vector
t = cellfun(@(x,y) x(y), trL, ind, 'unif',0);
trial_length = cell2mat(t);

t = [];
for itar = 1:length(targets)
    if isempty(temp_tr{itar})
        t = t;
    else
        t = cat(3,t,temp_tr{itar});
    end
end
trial_resp = t;

tt = [];
for itar = 1:length(targets)
   t = ones(1,size(temp_tr{itar},3))*targets(itar);
   tt = cat(2,tt,t);
end
trial_target = tt;



end



