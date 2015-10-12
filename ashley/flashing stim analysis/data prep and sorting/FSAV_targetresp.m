
%%
for i = 1:length(Dirs)
    ind = intersect(find(DirectionDeg == Dirs(i)),find(strcmp(trialOutcome,'success')));    
%     ind = intersect(find(DirectionDeg == Dirs(i)),find(strcmp(trialOutcome,'success') | strcmp(trialOutcome,'ignore')));    
%         ind = intersect(find(DirectionDeg == Dirs(i)),find(strcmp(trialOutcome,'ignore')));    
    if isempty(ind) == 1;
        dirTargetData{i} = [];
        dirTargetDataDF{i} = [];
        dirTargetDataDFoverF{i} = [];
        dirTargetInd{i} = ind;
%         dirTargetNorm{i} = [];
    elseif cTargetOn(1,ind(end))+40 > size(dataTC,1) 
        targetData = zeros(60,size(dataTC,2),length(ind));
        targetDataDF = zeros(60,size(dataTC,2),length(ind));
        targetDataDFoverF = zeros(60,size(dataTC,2),length(ind));
%         targetBL = zeros(length(Dirs),size(dataTC,2),length(ind));
%         targetNorm = zeros(60,size(dataTC,2),length(ind));
        for itrial = 1:length(ind)-1
            targetData(:,:,itrial) = dataTC(cTargetOn(ind(itrial))-19:cTargetOn(ind(itrial))+40,:);
            targetBaseline = dataTC(cLeverDown(ind(itrial))-29:cLeverDown(ind(itrial)),:);
            targetDataDF(:,:,itrial) = bsxfun(@minus, targetData(:,:,itrial), mean(targetBaseline,1));
            targetDataDFoverF(:,:,itrial) = bsxfun(@rdivide, targetDataDF(:,:,itrial), mean(targetBaseline,1));
%             targetBL(i,:,itrial) = mean(targetDataDFoverF(17:21,:,itrial),1);
%             targetNorm(:,:,itrial) = bsxfun(@minus,targetDataDFoverF(:,:,itrial),targetBL(i,:,itrial));
        end
        dirTargetData{i} = targetData;
        dirTargetDataDF{i} = targetDataDF;
        dirTargetDataDFoverF{i} = targetDataDFoverF;
        dirTargetInd{i} = ind;
%         dirTargetNorm{i} = targetNorm;        
    else
        targetData = zeros(60,size(dataTC,2),length(ind));
        targetDataDF = zeros(60,size(dataTC,2),length(ind));
        targetDataDFoverF = zeros(60,size(dataTC,2),length(ind));
%         targetBL = zeros(length(Dirs),size(dataTC,2),length(ind));
%         targetNorm = zeros(60,size(dataTC,2),length(ind));
        for itrial = 1:length(ind)
            targetData(:,:,itrial) = dataTC(cTargetOn(ind(itrial))-19:cTargetOn(ind(itrial))+40,:);
            targetBaseline = dataTC(cLeverDown(ind(itrial))-29:cLeverDown(ind(itrial)),:);
            targetDataDF(:,:,itrial) = bsxfun(@minus, targetData(:,:,itrial), mean(targetBaseline,1));
            targetDataDFoverF(:,:,itrial) = bsxfun(@rdivide, targetDataDF(:,:,itrial), mean(targetBaseline,1));
%             targetBL(i,:,itrial) = mean(targetDataDFoverF(17:21,:,itrial),1);
%             targetNorm(:,:,itrial) = bsxfun(@minus,targetDataDFoverF(:,:,itrial),targetBL(i,:,itrial));
        end
        dirTargetData{i} = targetData;
        dirTargetDataDF{i} = targetDataDF;
        dirTargetDataDFoverF{i} = targetDataDFoverF;
        dirTargetInd{i} = ind;
%         dirTargetNorm{i} = targetNorm;        
    end
    
end

%% mean response of cells to each target stim


