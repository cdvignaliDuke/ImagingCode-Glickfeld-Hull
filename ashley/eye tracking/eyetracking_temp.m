%% from raw data
eye = cat(3,Eye_data{:});
rad = sqrt(Area_temp./pi);
rad_norm = rad./nanmean(rad);
range_rad = [min(rad) max(rad)];

figure;
plot(rad_norm)
hold on
vline(find(isnan(rad_norm)),'r')

%find big and small pupils
[rad_sort rad_sort_ind] = sort(rad);
smallestpupils = rad_sort_ind(1:floor(0.001*length(rad)));

writetiff(eye(:,:,floor(28570/2)-50:floor(28570/2)+49),[fnout 'consecutivesmallpupil']);
writetiff(eye(:,:,floor(15999/2)-50:floor(15999/2)+49),[fnout 'smallestnon-consecutivepupil']);
writetiff(eye(:,:,floor(54066/2)-50:floor(54066/2)+49),[fnout 'biggestpupil']);

%find biggest change in pupil size
[rad_sort rad_sort_ind] = sort(abs(diff(rad)));
rad_sort_ind = rad_sort_ind(1:length(rad_sort_ind)-sum(isnan(rad_sort))-1);
biggestchange_ind = rad_sort_ind(end-sum(isnan(rad_sort))-1)-5:rad_sort_ind(end-sum(isnan(rad_sort))-1)+4;
writetiff(eye(:,:,floor(biggestchange_ind/2)),[fnout 'biggestpupilchange'])
%% time-courses
figure;
subplot(1,2,1)
plot(rad_mat_down_base(:,b1Ix))
subplot(1,2,2)
plot(rad_mat_down_base(:,b2Ix))

rad_norm_v = nanmean(rad_mat_down_base(:,b1Ix),2);
rad_norm_a = nanmean(rad_mat_down_base(:,b2Ix),2);

figure;
plot(rad_norm_v,'g')
hold on
plot(rad_norm_a,'k')

%% baseline variability?

figure;
suptitle('baseline raw value')
subplot(1,3,1)
[N_v,X_v] = hist(nanmean(rad_mat_down(1:30,b1Ix),1));
hist_v = bar(X_v,N_v);
hist_v.EdgeColor = [1,1,1];
hist_v.FaceColor = 'g';
hist_v.FaceAlpha = .5;
hold on
[N_a,X_a] = hist(nanmean(rad_mat_down(1:30,b2Ix),1));
hist_a = bar(X_a,N_a);
hist_a.EdgeColor = [1,1,1];
hist_a.FaceColor = 'k';
hist_a.FaceAlpha = .5;

subplot(1,3,2)
[N_v,X_v] = hist(nanmean(rad_mat_down(1:15,b1Ix),1));
hist_v = bar(X_v,N_v);
hist_v.EdgeColor = [1,1,1];
hist_v.FaceColor = 'g';
hist_v.FaceAlpha = .5;
hold on
[N_a,X_a] = hist(nanmean(rad_mat_down(1:15,b2Ix),1));
hist_a = bar(X_a,N_a);
hist_a.EdgeColor = [1,1,1];
hist_a.FaceColor = 'k';
hist_a.FaceAlpha = .5;

subplot(1,3,3)
[N_v,X_v] = hist(nanmean(rad_mat_down(16:30,b1Ix),1));
hist_v = bar(X_v,N_v);
hist_v.EdgeColor = [1,1,1];
hist_v.FaceColor = 'g';
hist_v.FaceAlpha = .5;
hold on
[N_a,X_a] = hist(nanmean(rad_mat_down(16:30,b2Ix),1));
hist_a = bar(X_a,N_a);
hist_a.EdgeColor = [1,1,1];
hist_a.FaceColor = 'k';
hist_a.FaceAlpha = .5;

% ttest start and end of baseline period
early_base = rad_mat_down(1:5,:);
late_base = rad_mat_down(end-4:end,:);

[H,P] = ttest(early_base,late_base);

% plot stable and unstable baseline trials.
figure;
subplot(3,2,1)
plot(rad_mat_down_base(:,find(H)));
title('baseline varies')
subplot(3,2,2)
plot(rad_mat_down_base(:,find(H==0)))
title('baseline flat')
subplot(3,2,3)
plot(nanmean(rad_mat_down(:,intersect(find(H),b1Ix)),2),'g')
hold on
plot(nanmean(rad_mat_down(:,intersect(find(H),b2Ix)),2),'k')
subplot(3,2,4)
plot(nanmean(rad_mat_down(:,intersect(find(H==0),b1Ix)),2),'g')
hold on
plot(nanmean(rad_mat_down(:,intersect(find(H==0),b2Ix)),2),'k')
subplot(3,2,5)
plot(nanmean(rad_mat_down_base(:,intersect(find(H),b1Ix)),2),'g')
hold on
plot(nanmean(rad_mat_down_base(:,intersect(find(H),b2Ix)),2),'k')
subplot(3,2,6)
plot(nanmean(rad_mat_down_base(:,intersect(find(H==0),b1Ix)),2),'g')
hold on
plot(nanmean(rad_mat_down_base(:,intersect(find(H==0),b2Ix)),2),'k')



%% quantification
vistrials = 1;
audtrials = 2;
alltrialoutcome = 1;

pressalign = figure;
for iplot = 1:6
    figure(pressalign)
    subplot(3,2,iplot)
    hold on
    plot([-10:10],[-10:10],'k--')
    hold on
    vline(1,'k');
    hline(1,'k');
    xlim([0 2]);
    ylim([0 2]);
    xlabel('visual');
    ylabel('auditory');
    axis square
end
targetalign = figure;

for ialign = 1:3;
for imouse = 1:nMice
    for iexp = 1:size(mouse.expt,2)
        
        %plot absolute radius value
        radius_trans_v =  mouse(imouse).expt(iexp).align(ialign).av(vistrials).outcome(alltrialoutcome).avg_rad_trans; 
        radius_trans_a =  mouse(imouse).expt(iexp).align(ialign).av(audtrials).outcome(alltrialoutcome).avg_rad_trans; 
        figure(pressalign)
        subplot(3,2,1)
        hold on
        rad_trans_plot = errorbarxy(radius_trans_v(1),radius_trans_a(1),radius_trans_v(2),radius_trans_a(2),{'bo','b','b'});
        
        %plot normalized radius value
        radius_trans_v =  mouse(imouse).expt(iexp).align(ialign).av(vistrials).outcome(alltrialoutcome).avg_rad_trans ./ mouse(imouse).expt(iexp).align(ialign).av(vistrials).outcome(alltrialoutcome).avg_rad_pre;
        radius_trans_a =  mouse(imouse).expt(iexp).align(ialign).av(audtrials).outcome(alltrialoutcome).avg_rad_trans ./ mouse(imouse).expt(iexp).align(ialign).av(audtrials).outcome(alltrialoutcome).avg_rad_pre;
        figure(pressalign)
        subplot(3,2,2)
        hold on
        rad_trans_plot = scatter(radius_trans_v(1),radius_trans_a(1));
        
        %absolute horizontal pos
        horz_trans_v = mouse(imouse).expt(iexp).align(ialign).av(vistrials).outcome(alltrialoutcome).avg_hor_trans; 
        horz_trans_a = mouse(imouse).expt(iexp).align(ialign).av(audtrials).outcome(alltrialoutcome).avg_hor_trans; 
        figure(pressalign)
        subplot(3,2,3)
        hold on
        rad_trans_plot = errorbarxy(horz_trans_v(1),horz_trans_a(1),horz_trans_v(2),horz_trans_a(2),{'bo','b','b'});
        
        
        
        