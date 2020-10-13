ntrials = length(celleqel2mat_padded(input.tGratingContrast));
data_avg = squeeze(mean(mean(data,1),2));
data_avg_top = squeeze(mean(mean(data(1:130,:,:),1),2));
data_avg_bottom = squeeze(mean(mean(data(130:end,:,:),1),2));
% nOn = input.nScansOn;
% nOff = input.nScansOff;
% cFirstStim = celleqel2mat_padded(input.cFirstStim);
% cTarget = celleqel2mat_padded(input.cTargetOn);
% cStimOn = celleqel2mat_padded(input.cStimOn);
cStimOff = celleqel2mat_padded(input.cStimOff);
data_mat = zeros(20,ntrials);
data_mat_top = zeros(20,ntrials);
data_mat_bottom = zeros(20,ntrials);
% start = 1;
% for itrial = 1:ntrials
%     cOn = start+nOff;
%     data_mat(:,itrial) = data_avg(cOn-10:cOn+9,:);
%     data_mat_top(:,itrial) = data_avg_top(cOn-10:cOn+9,:);
%     data_mat_bottom(:,itrial) = data_avg_bottom(cOn-10:cOn+9,:);
%     start = nOn+cOn;
% end

% for itrial = 1:ntrials
%     data_mat(:,itrial) = data_avg(cFirstStim(itrial)-10:cFirstStim(itrial)+9,:);
%     data_mat_top(:,itrial) = data_avg_top(cFirstStim(itrial)-10:cFirstStim(itrial)+9,:);
%     data_mat_bottom(:,itrial) = data_avg_bottom(cFirstStim(itrial)-10:cFirstStim(itrial)+9,:);
% end

% for itrial = 1:ntrials
%     data_mat(:,itrial) = data_avg(cTarget(itrial)-10:cTarget(itrial)+9,:);
%     data_mat_top(:,itrial) = data_avg_top(cTarget(itrial)-10:cTarget(itrial)+9,:);
%     data_mat_bottom(:,itrial) = data_avg_bottom(cTarget(itrial)-10:cTarget(itrial)+9,:);
% end

% for itrial = 1:ntrials
%     data_mat(:,itrial) = data_avg(cStimOn(itrial)-10:cStimOn(itrial)+9,:);
%     data_mat_top(:,itrial) = data_avg_top(cStimOn(itrial)-10:cStimOn(itrial)+9,:);
%     data_mat_bottom(:,itrial) = data_avg_bottom(cStimOn(itrial)-10:cStimOn(itrial)+9,:);
% end


for itrial = 1:ntrials
    data_mat(:,itrial) = data_avg(cStimOff(itrial)-10:cStimOff(itrial)+9,:);
    data_mat_top(:,itrial) = data_avg_top(cStimOff(itrial)-10:cStimOff(itrial)+9,:);
    data_mat_bottom(:,itrial) = data_avg_bottom(cStimOff(itrial)-10:cStimOff(itrial)+9,:);
end

data_mat_sub = data_mat - mean(data_mat(1:10,:),1);
data_mat_sub_top = data_mat_top - mean(data_mat_top(1:10,:),1);
data_mat_sub_bottom = data_mat_bottom - mean(data_mat_bottom(1:10,:),1);
tt = [-10:9] .*(1000/30);

% figure;
% suptitle('170710- VisStimRet')
% subplot(2,2,1)
% plot(tt(8:14),data_mat_sub(8:14,:))
% vline(0)
% subplot(2,2,2)
% plot(tt,mean(data_mat_sub,2))
% title('Whole FOV')
% vline(0)
% subplot(2,2,3)
% plot(tt,mean(data_mat_sub_top,2))
% vline(0)
% title('Top')
% subplot(2,2,4)
% plot(tt,mean(data_mat_sub_bottom,2))
% vline(0)
% title('Bottom')

figure;
subplot(2,2,1)
plot(tt(8:14),data_mat_sub(8:14,:))
vline(0)
subplot(2,2,2)
plot(tt,mean(data_mat_sub,2))
title('Whole FOV')
vline(0)
subplot(2,2,3)
plot(tt,mean(data_mat_sub_top,2))
vline(0)
title('Top')
subplot(2,2,4)
plot(tt,mean(data_mat_sub_bottom,2))
vline(0)
title('Bottom')
suptitle('180307- Lego2PFrames- StimOff')