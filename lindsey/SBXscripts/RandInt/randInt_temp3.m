data_dfof_trialavg = squeeze(mean(bsxfun(@minus, mean(data_dfof(29:34,:,:,:),1),mean(data_dfof(19:24,:,:,:),1)),4));
data_dfof_trialavg_thresh = data_dfof_trialavg;
data_dfof_trialavg_thresh(data_dfof_trialavg <0) = 0;
data_dfof_pp = data_dfof_trialavg_thresh(:,2)./data_dfof_trialavg_thresh(:,1);
data_dfof_ss = mean(data_dfof_trialavg_thresh(:,3:5),2)./data_dfof_trialavg_thresh(:,1);
figure;
subplot(2,2,1)
hist(data_dfof_pp(good_ind))
title('Paired pulse')
xlim([0 2])
subplot(2,2,2)
hist(data_dfof_ss(good_ind))
title('Steady state')
xlim([0 2])
subplot(2,2,3)
scatter(data_dfof_trialavg_thresh(good_ind,1),data_dfof_pp(good_ind), 'ok')
xlabel('Base amp')
ylabel('Paired pulse')
ylim([0 2])
subplot(2,2,4)
scatter(data_dfof_trialavg_thresh(good_ind,1),data_dfof_ss(good_ind), 'ok')
xlabel('Base amp')
ylabel('Steady state')
ylim([0 2])

tt = [-20:19].*input.frameRateHz;
ind = find(baseDir == dirs(1));
[a b] = ttest(squeeze(mean(data_dfof_trial(29:34,:,1,ind),1))', squeeze(mean(data_dfof_trial(19:24,:,1,ind),1))', 'tail', 'right');
good_base_ind = find(a);
[n n2] = subplotn(nCells);
figure;
for iC = 1:nCells
    subplot(n,n2,iC)
    plot(tt, mean(data_dfof_trial(:,iC,1,:),4))
    hold on
    for ioff = 1:noff
        ind = find(tFramesOff(:,1) == offs(ioff));
        plot(tt, mean(data_dfof_trial(:,iC,2,ind),4))
    end
    if find(good_ind == iC)
        title(num2str(dirs(max_dir(iC))))
    end
end

nCells = size(data_dfof,2);
data_dfof_pp = zeros(nCells,noff);
data_dfof_ss = zeros(nCells,noff);
data_dfof_trialavg = zeros(nCells,nCyc,noff);
for ioff = 1:noff
    figure;
    data_dfof_trialavg(:,1,ioff) = squeeze(mean(bsxfun(@minus, mean(data_dfof(29:34,:,1,:),1),mean(data_dfof(19:24,:,1,:),1)),4));
    for icyc = 2:nCyc-1
        ind = find(tFramesOff(:,icyc-1) == offs(ioff));
        data_dfof_trialavg(:,icyc,ioff) = squeeze(mean(bsxfun(@minus, mean(data_dfof(29:34,:,icyc,ind),1),mean(data_dfof(19:24,:,icyc,ind),1)),4));
    end
    data_dfof_trialavg_thresh = data_dfof_trialavg;
    data_dfof_trialavg_thresh(data_dfof_trialavg <0) = 0;
    data_dfof_pp(:,ioff) = data_dfof_trialavg_thresh(:,2,ioff)./data_dfof_trialavg_thresh(:,1,ioff);
    data_dfof_ss(:,ioff) = mean(data_dfof_trialavg_thresh(:,3:5,ioff),2)./data_dfof_trialavg_thresh(:,1,ioff);
    subplot(2,2,1)
    hist(data_dfof_pp(good_ind,ioff),[0:0.2:2.6])
    title(['Paired pulse- mean: ' num2str(chop(mean(data_dfof_pp(good_ind,ioff),1),2)) ])
    xlim([0 2.5])
    subplot(2,2,2)
    hist(data_dfof_ss(good_ind,ioff),[0:0.2:2.6])
    title(['Steady state- mean: ' num2str(chop(mean(data_dfof_ss(good_ind,ioff),1),2)) ])
    xlim([0 2.5])
    subplot(2,2,3)
    scatter(data_dfof_trialavg_thresh(good_ind,1),data_dfof_pp(good_ind,ioff), 'ok')
    xlabel('Base amp')
    ylabel('Paired pulse')
    ylim([0 2])
    subplot(2,2,4)
    scatter(data_dfof_trialavg_thresh(good_ind,1),data_dfof_ss(good_ind,ioff), 'ok')
    xlabel('Base amp')
    ylabel('Steady state')
    ylim([0 2])
    suptitle([num2str(offs(ioff)*input.frameRateHz) 'ms'])
end

figure; 
for i = 1:noff
    subplot(2,3,i)
    scatter(data_dfof_trialavg_thresh(good_ind,2,i),data_dfof_trialavg_thresh(good_ind,1,i),'o')
    hold on
    plot(0:.1:.3,0:.1:.3,'-')
    xlim([0 0.3])
    ylim([0 0.3])
    axis square
    title([num2str(offs(i)*input.frameRateHz) 'ms'])
    xlabel('Pulse 2')
    ylabel('Pulse 1')
    subplot(2,3,i+noff)
    scatter(mean(data_dfof_trialavg_thresh(good_ind,3:5,i),2),data_dfof_trialavg_thresh(good_ind,1,i),'o')
    hold on
    plot(0:.1:.3,0:.1:.3,'-')
    xlim([0 0.3])
    ylim([0 0.3])
    axis square
    title([num2str(offs(i)*input.frameRateHz) 'ms'])
    xlabel('Pulse 3-5')
    ylabel('Pulse 1')
end

figure;
subplot(2,2,1)
scatter(data_dfof_pp(good_ind,1),data_dfof_pp(good_ind,3))
xlabel([num2str(offs(1)*input.frameRateHz) 'ms'])
ylabel([num2str(offs(3)*input.frameRateHz) 'ms'])
title('Paired pulse')
hold on
plot(0:.1:5,0:.1:5)
xlim([0 5])
ylim([0 5])
axis square
subplot(2,2,2)
scatter(data_dfof_ss(good_ind,1),data_dfof_ss(good_ind,3))
xlabel([num2str(offs(1)*input.frameRateHz) 'ms'])
ylabel([num2str(offs(3)*input.frameRateHz) 'ms'])
title('Steady state')
hold on
plot(0:.1:2,0:.1:2)
xlim([0 2])
ylim([0 2])
axis square

