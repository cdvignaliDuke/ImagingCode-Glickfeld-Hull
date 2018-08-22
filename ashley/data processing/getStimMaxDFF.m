function stimMaxDFF = getStimMaxDFF(data,trialStart,trialType,basewin,respwin);

types = unique(trialType);
nTypes = length(types);
[ypix,xpix,nfr] = size(data);
lastFrame = respwin(end);
nBaseFrame = basewin(end);
nStimFrame = lastFrame - nBaseFrame;

if trialStart(end)+(lastFrame-nBaseFrame) > nfr
    ind = find(trialStart+(lastFrame-nBaseFrame) > nfr,1);
    trialStart = trialStart(1:(ind-1));
    trialType = trialType(1:(ind-1));
elseif length(trialStart) ~= length(trialType)
    if length(trialStart) > length(trialType)
        trialStart = trialStart(1:length(trialType));
    else
        trialType = trialType(1:length(trialStart));     
    end
end

F = cell(1,nTypes);
for i = 1:nTypes
    ind = find(trialType == i);
    d = nan(ypix,xpix,length(ind),lastFrame);
    for itrial = 1:length(ind)
        trialInd = trialStart(ind(itrial));
        d(:,:,itrial,:) = data(:,:,(trialInd-nBaseFrame+1):(trialInd+nStimFrame));
    end
    F{i} = d;
end

F0 = cellfun(@(x) mean(x(:,:,:,basewin),4),F,'unif',0);
dFF = cellfun(@(x,y) (x-y)./y,F,F0,'unif',0);

respDFF = cellfun(@(x) mean(x(:,:,:,respwin),4), dFF,'unif',0);
maxDFF = nan(ypix,xpix,nTypes);
for i = 1:nTypes
    maxDFF(:,:,i) = max(respDFF{i},[],3);
end

stimMaxDFF = maxDFF;

end