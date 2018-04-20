Allmice = {'527','549','559','563'};
num_subj = size(Allmice,2);
Thresh = NaN(num_subj,3);
AbsThresh = NaN(num_subj,3);
FA = NaN(num_subj,3);
Thresh_off = NaN(num_subj,3,3);
AbsThresh_off = NaN(num_subj,3,3);
FA_off = NaN(num_subj,3,3);

FA_first = NaN(num_subj,3);
Thresh_first = NaN(num_subj,3);
AbsThresh_first = NaN(num_subj,3);
FA_last = NaN(num_subj,3);
Thresh_last = NaN(num_subj,3);
AbsThresh_last = NaN(num_subj,3);




for i_mice = 1:num_subj
    name = Allmice{i_mice};
    out = load(name);
    Thresh(i_mice,:) = out.Threshall;
    FA(i_mice,:) = out.RealFA_B;
    AbsThresh(i_mice,:) = out.AbsThreshall;
    
    Thresh_first(i_mice,:) = out.first.Threshall;
    Thresh_last(i_mice,:) = out.last.Threshall;
    AbsThresh_first(i_mice,:) = out.first.AbsThreshall;
    AbsThresh_last(i_mice,:) = out.last.AbsThreshall;
    FA_first(i_mice,:) = [out.first.FA{1}.all.FA out.first.FA{2}.all.FA out.first.FA{3}.all.FA];
    FA_last(i_mice,:) = [out.last.FA{1}.all.FA out.last.FA{2}.all.FA out.last.FA{3}.all.FA];
    
    Thresh_off(i_mice,:,:) = out.Thresh;
    AbsThresh_off(i_mice,:,:) = out.AbsThresh;
    for i = 1:3
    FA_off(i_mice,i,:) = (out.output.FA{i}.FA)';
    end 
end 
B_orien = out.B_ori;
Off = out.Off;
%% plot the threshold and FA rate difference between different base orientation
figure
subplot(2,2,1)
NormThresh = Thresh./repmat(Thresh(:,1),1,3);
Mean_norm = mean(NormThresh,1);
sem_norm = std(NormThresh,[],1)./sqrt(num_subj);
h=errorbar([2 1 3],Mean_norm,sem_norm,'k');
h.LineStyle = 'none';
hold on 
scatter([2 1 3],Mean_norm,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'sizeData',50 )

ylim([0 2])
xlim([0.5 3.5])
hline(1,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'15','0','165'},'TickDir','out')
xlabel('Base Orien (Deg)')
ylabel('Thresh delta orien (deg)')

subplot(2,2,2)
NormThresh = AbsThresh./repmat(AbsThresh(:,1),1,3);
Mean_norm = mean(NormThresh,1);
sem_norm = std(NormThresh,[],1)./sqrt(num_subj);
h=errorbar([2 1 3],Mean_norm,sem_norm,'k');
h.LineStyle = 'none';
hold on 
scatter([2 1 3],Mean_norm,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'sizeData',50 )

ylim([0 2])
xlim([0.5 3.5])
hline(1,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'15','0','165'},'TickDir','out')
xlabel('Base Orien (Deg)')
ylabel('Thresh abs orien (deg)')

subplot(2,2,3)
NormThresh = FA./repmat(FA(:,1),1,3);
Mean_norm = mean(NormThresh,1);
sem_norm = std(NormThresh,[],1)./sqrt(num_subj);
h=errorbar([2 1 3],Mean_norm,sem_norm,'k');
h.LineStyle = 'none';
hold on 
scatter([2 1 3],Mean_norm,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'sizeData',50 )

ylim([0 2.5])
xlim([0.5 3.5])
hline(1,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'15','0','165'},'TickDir','out')
xlabel('Base Orien (Deg)')
ylabel('FA rate')


%% plot the threshold and FA rate difference between different base orientation, overlay different animal's performance
colormice = lines(4);
figure
subplot(2,2,1)
NormThresh = Thresh - repmat(Thresh(:,1),1,3);
Mean_norm = mean(NormThresh,1);
sem_norm = std(NormThresh,[],1)./sqrt(num_subj);
for i = 1:num_subj
    scatter([0.8 3.2],NormThresh(i,[2 3]),'MarkerEdgeColor',colormice(i,:),'MarkerFaceColor',colormice(i,:),'sizeData',15 )
    hold on
end 
h=errorbar([2 1 3],Mean_norm,sem_norm,'k');
h.LineStyle = 'none';
hold on 
scatter([2 1 3],Mean_norm,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'sizeData',50 )
text(0.7,18,'Orien Change')
ylim([-20 20])
xlim([0.5 3.5])
hline(0,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'15','0','165'},'TickDir','out')
xlabel('Base Orien (Deg)')
ylabel('\Delta Threshold (deg)')
% statistical test
stat_temp = zeros(4,3);
stat_temp(:,1) = Thresh(:,2);
stat_temp(:,2) = Thresh(:,1);
stat_temp(:,3) = Thresh(:,3);
[p,tbl,stats] = anova1(stat_temp);
[results,means] = multcompare(stats);
%

subplot(2,2,2)
NormThresh = AbsThresh - repmat(AbsThresh(:,1),1,3);
Mean_norm = mean(NormThresh,1);
sem_norm = std(NormThresh,[],1)./sqrt(num_subj);
for i = 1:num_subj
    scatter([0.8 3.2],NormThresh(i,[2 3]),'MarkerEdgeColor',colormice(i,:),'MarkerFaceColor',colormice(i,:),'sizeData',15)
    hold on
end 

h=errorbar([2 1 3],Mean_norm,sem_norm,'k');
h.LineStyle = 'none';
hold on 
scatter([2 1 3],Mean_norm,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'sizeData',50 )
text(0.7,18,'Absolute Orien')
ylim([-20 20])
xlim([0.5 3.5])
hline(0,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'15','0','165'},'TickDir','out')
xlabel('Base Orien (Deg)')
ylabel('\Delta Threshold (deg)')

subplot(2,2,3)
NormThresh = FA - repmat(FA(:,1),1,3);
Mean_norm = mean(NormThresh,1);
sem_norm = std(NormThresh,[],1)./sqrt(num_subj);
for i = 1:num_subj
    scatter([0.8 3.2],NormThresh(i,[2 3]),'MarkerEdgeColor',colormice(i,:),'MarkerFaceColor',colormice(i,:),'sizeData',15)
    hold on
end 
h=errorbar([2 1 3],Mean_norm,sem_norm,'k');
h.LineStyle = 'none';
hold on 
scatter([2 1 3],Mean_norm,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'sizeData',50 )

 ylim([-0.05 0.1])
xlim([0.5 3.5])
hline(0,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'15','0','165'},'TickDir','out')
xlabel('Base Orien (Deg)')
ylabel('\Delta FA')

%% plot threshold 
figure
NormThresh = Thresh - repmat(Thresh(:,1),1,3);
Mean_norm = mean(NormThresh,1);
sem_norm = std(NormThresh,[],1)./sqrt(num_subj);
for i = 1:num_subj
    scatter([1 3],NormThresh(i,[2 3]),'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'sizeData',10)
    hold on
end 
h=errorbar([2 1 3],Mean_norm,sem_norm,'k');
h.LineStyle = 'none';
hold on 
scatter([2 1 3],Mean_norm,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'sizeData',50 )

ylim([-20 20])
xlim([0.7 3.3])
hline([-15 15],'k:')
hline(0,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'15','0','165'},'TickDir','out')
xlabel('Base Orien (Deg)')
ylabel('\Delta Threshold (deg)')
%% plot the off dependence of threshold effects

colors(1,:) = [0 0 0];
colors(2,:) = [1 0.2 0]; % red 15 deg
colors(3,:) = [0.2 0.6 0]; % green 165 deg
subplot(3,3,1)

temp = Thresh_off - repmat(Thresh_off(:,1,:),1,3,1);
mean_temp =squeeze( mean(temp,1));
sem_temp = squeeze(std(temp,[],1))./sqrt(num_subj);
for i = 1:3
    h = errorbar([1 2 3],mean_temp(i,:),sem_temp(i,:),'color',colors(i,:));
    h.LineStyle = 'none';
    hold on
    scatter([1 2 3],mean_temp(i,:),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1 1 1],'sizeData',50)
    
end 
text(0.7, 18, 'Orien Change')
ylim([-20 20])
xlim([0.5 3.5])
hline(0,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('\Delta Threshold(deg)')


subplot(3,3,2)

temp = Thresh_off./repmat(Thresh_off(:,:,1),1,1,3);
mean_temp =squeeze( mean(temp,1));
sem_temp = squeeze(std(temp,[],1))./sqrt(num_subj);
for i = 1:3
    h = errorbar([1 2 3],mean_temp(i,:),sem_temp(i,:),'color',colors(i,:));
    h.LineStyle = 'none';
    hold on
    scatter([1 2 3],mean_temp(i,:),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1 1 1],'sizeData',50)
    text(0.7,1-0.2*i,[num2str(B_orien(i)) 'Deg'],'color',colors(i,:))
end 

ylim([0 1.1])
xlim([0.5 3.5])
hline(1,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Norm. Threshold')

subplot(3,3,3)

temp = Thresh_off./repmat(Thresh_off(:,1,1),1,3,3);
mean_temp =squeeze( mean(temp,1));
sem_temp = squeeze(std(temp,[],1))./sqrt(num_subj);
for i = 1:3
    h = errorbar([1 2 3],mean_temp(i,:),sem_temp(i,:),'color',colors(i,:));
    h.LineStyle = 'none';
    hold on
    scatter([1 2 3],mean_temp(i,:),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1 1 1],'sizeData',50)
    
end 
text(0.7, 1.8, 'Orien Change')
ylim([0 2])
xlim([0.5 3.5])
hline(1,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Norm. Threshold')

subplot(3,3,4)

temp = AbsThresh_off - repmat(AbsThresh_off(:,1,:),1,3,1);
mean_temp =squeeze( mean(temp,1));
sem_temp = squeeze(std(temp,[],1))./sqrt(num_subj);
for i = 1:3
    h = errorbar([1 2 3],mean_temp(i,:),sem_temp(i,:),'color',colors(i,:));
    h.LineStyle = 'none';
    hold on
    scatter([1 2 3],mean_temp(i,:),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1 1 1],'sizeData',50)
    
end 
text(0.7, 18, 'Absolute Orien')
ylim([-20 20])
xlim([0.5 3.5])
hline(0,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('\Delta Threshold(deg)')

subplot(3,3,5)

temp = AbsThresh_off./repmat(AbsThresh_off(:,:,1),1,1,3);
mean_temp =squeeze( mean(temp,1));
sem_temp = squeeze(std(temp,[],1))./sqrt(num_subj);
for i = 1:3
    h = errorbar([1 2 3],mean_temp(i,:),sem_temp(i,:),'color',colors(i,:));
    h.LineStyle = 'none';
    hold on
    scatter([1 2 3],mean_temp(i,:),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1 1 1],'sizeData',50)
    
end 

ylim([0 1.1])
xlim([0.5 3.5])
hline(1,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Norm. Threshold')

subplot(3,3,6)

temp = AbsThresh_off./repmat(AbsThresh_off(:,1,1),1,3,3);
mean_temp =squeeze( mean(temp,1));
sem_temp = squeeze(std(temp,[],1))./sqrt(num_subj);
for i = 1:3
    h = errorbar([1 2 3],mean_temp(i,:),sem_temp(i,:),'color',colors(i,:));
    h.LineStyle = 'none';
    hold on
    scatter([1 2 3],mean_temp(i,:),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1 1 1],'sizeData',50)
    
end 
text(0.7, 1.8, 'Absolute Orien')
ylim([0 2])
xlim([0.5 3.5])
hline(1,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Norm. Threshold')

subplot(3,3,7)

temp = FA_off - repmat(FA_off(:,1,:),1,3,1);
mean_temp =squeeze( mean(temp,1));
sem_temp = squeeze(std(temp,[],1))./sqrt(num_subj);
for i = 1:3
    h = errorbar([1 2 3],mean_temp(i,:),sem_temp(i,:),'color',colors(i,:));
    h.LineStyle = 'none';
    hold on
    scatter([1 2 3],mean_temp(i,:),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1 1 1],'sizeData',50)
    
end 

ylim([-0.04 0.08])
xlim([0.5 3.5])
hline(0,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('\Delta FA')

subplot(3,3,8)

temp = FA_off./repmat(FA_off(:,:,1),1,1,3);
mean_temp =squeeze( mean(temp,1));
sem_temp = squeeze(std(temp,[],1))./sqrt(num_subj);
for i = 1:3
    h = errorbar([1 2 3],mean_temp(i,:),sem_temp(i,:),'color',colors(i,:));
    h.LineStyle = 'none';
    hold on
    scatter([1 2 3],mean_temp(i,:),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1 1 1],'sizeData',50)
   
end 

ylim([0 4])
xlim([0.5 3.5])
hline(1,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Norm. FA')

subplot(3,3,9)

temp = FA_off./repmat(FA_off(:,1,1),1,3,3);
mean_temp =squeeze( mean(temp,1));
sem_temp = squeeze(std(temp,[],1))./sqrt(num_subj);
for i = 1:3
    h = errorbar([1 2 3],mean_temp(i,:),sem_temp(i,:),'color',colors(i,:));
    h.LineStyle = 'none';
    hold on
    scatter([1 2 3],mean_temp(i,:),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1 1 1],'sizeData',50)
   
end 

ylim([0 6])
xlim([0.5 3.5])
hline(1,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Norm. FA')

%% plot the FA rate

figure
NormThresh = FA - repmat(FA(:,1),1,3);
Mean_norm = mean(NormThresh,1);
sem_norm = std(NormThresh,[],1)./sqrt(num_subj);
for i = 1:num_subj
    scatter([1 3],NormThresh(i,[2 3]),'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'sizeData',10)
    hold on
end 
h=errorbar([2 1 3],Mean_norm,sem_norm,'k');
h.LineStyle = 'none';
hold on 
scatter([2 1 3],Mean_norm,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'sizeData',50 )

 ylim([-0.02 0.08])
xlim([0.7 3.3])

hline(0,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'15','0','165'},'TickDir','out')
xlabel('Base Orien (Deg)')
ylabel('\Delta FA rate')
%% summary the norm threshold to 0 deg 250 condition
figure
temp = Thresh_off./repmat(Thresh_off(:,1,1),1,3,3);
mean_temp =squeeze( mean(temp,1));
sem_temp = squeeze(std(temp,[],1))./sqrt(num_subj);
for i = 1:num_subj
    for i_base = 1:3
    for i_off = 1:3
    scatter(i_off,squeeze(temp(i,i_base,i_off)),'MarkerEdgeColor',colors(i_base,:),'MarkerFaceColor',colors(i_base,:),'sizeData',10)
    hold on
    end 
    end 
end

for i = 1:3
    h = errorbar([1 2 3],mean_temp(i,:),sem_temp(i,:),'color',colors(i,:));
    h.LineStyle = 'none';
    hold on
    scatter([1 2 3],mean_temp(i,:),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1 1 1],'sizeData',50)
    
end 

%  ylim([0 2])
xlim([0.5 3.5])
hline(1,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Norm. Threshold')

%% summary the norm FA rate to 0 deg 250 condition
figure
temp = FA_off./repmat(FA_off(:,1,1),1,3,3);
mean_temp =squeeze( mean(temp,1));
sem_temp = squeeze(std(temp,[],1))./sqrt(num_subj);
for i = 1:num_subj
    for i_base = 1:3
    for i_off = 1:3
    scatter(i_off,squeeze(temp(i,i_base,i_off)),'MarkerEdgeColor',colors(i_base,:),'MarkerFaceColor',colors(i_base,:),'sizeData',10)
    hold on
    end 
    end 
end

for i = 1:3
    h = errorbar([1 2 3],mean_temp(i,:),sem_temp(i,:),'color',colors(i,:));
    h.LineStyle = 'none';
    hold on
    scatter([1 2 3],mean_temp(i,:),'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1 1 1],'sizeData',50)
    
end 

%  ylim([0 2])
xlim([0.5 3.5])
hline(1,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Norm. Threshold')
%% summarize the training stuff first third vs last third
figure
subplot(1,2,1)
NormThresh = Thresh_first - repmat(Thresh_first(:,1),1,3);
Mean_norm = mean(NormThresh,1);
sem_norm = std(NormThresh,[],1)./sqrt(num_subj);
for i = 1:num_subj
    scatter([1 3],NormThresh(i,[2 3]),'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'sizeData',10)
    hold on
end 
h=errorbar([2 1 3],Mean_norm,sem_norm,'k');
h.LineStyle = 'none';
hold on 
scatter([2 1 3],Mean_norm,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'sizeData',50 )

 ylim([-20 20])
xlim([0.7 3.3])
hline([-15 15],'k:')
hline(0,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'15','0','165'},'TickDir','out')
xlabel('Base Orien (Deg)')
ylabel('\Delta Threshold (deg)')
subplot(1,2,2)
NormThresh = Thresh_last - repmat(Thresh_last(:,1),1,3);
Mean_norm = mean(NormThresh,1);
sem_norm = std(NormThresh,[],1)./sqrt(num_subj);
for i = 1:num_subj
    scatter([1 3],NormThresh(i,[2 3]),'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'sizeData',10)
    hold on
end 
h=errorbar([2 1 3],Mean_norm,sem_norm,'k');
h.LineStyle = 'none';
hold on 
scatter([2 1 3],Mean_norm,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'sizeData',50 )

 ylim([-20 20])
xlim([0.7 3.3])
hline([-15 15],'k:')
hline(0,'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'15','0','165'},'TickDir','out')
xlabel('Base Orien (Deg)')
ylabel('\Delta Threshold (deg)')

%% plot the first third and last third in the same axis 
NormThresh_first = Thresh_first - repmat(Thresh_first(:,1),1,3);
NormThresh_last = Thresh_last - repmat(Thresh_last(:,1),1,3);
mean_first = mean(NormThresh_first,1);
sem_first = std(NormThresh_first,[],1)./sqrt(num_subj);
mean_last = mean(NormThresh_last,1);
sem_last = std(NormThresh_last,[],1)./sqrt(num_subj);
figure
subplot(1,2,1)
for i = 1:num_subj
    plot([0.9 1.1],[NormThresh_first(i,2) NormThresh_last(i,2)], 'color',colormice(i,:),'Marker','o')
    hold on
    
end 

h=errorbar([0.9 1.1],[mean_first(2) mean_last(2)],[sem_first(2) sem_last(2)],'k','Marker','o');
h.LineStyle = 'none';

for i = 1:num_subj
    plot([1.9 2.1],[NormThresh_first(i,3) NormThresh_last(i,3)], 'color',colormice(i,:),'Marker','o')
    hold on
    
end 

h=errorbar([1.9 2.1],[mean_first(3) mean_last(3)],[sem_first(3) sem_last(3)],'k','Marker','o');
h.LineStyle = 'none';

ylim([-20 20])
xlim([0.7 2.3])
hline([-15 15],'k:')
hline(0,'k:')
set(gca,'XTick',[1 2],'XTickLabel',{'15','-15'},'TickDir','out')
xlabel('Base Orien (Deg)')
ylabel('\Delta Threshold (deg)')

% plot the FA rate
NormThresh_first = FA_first - repmat(FA_first(:,1),1,3);
NormThresh_last = FA_last - repmat(FA_last(:,1),1,3);
mean_first = mean(NormThresh_first,1);
sem_first = std(NormThresh_first,[],1)./sqrt(num_subj);
mean_last = mean(NormThresh_last,1);
sem_last = std(NormThresh_last,[],1)./sqrt(num_subj);

subplot(1,2,2)
for i = 1:num_subj
    plot([0.9 1.1],[NormThresh_first(i,2) NormThresh_last(i,2)], 'color',colormice(i,:),'Marker','o')
    hold on
    
end 

h=errorbar([0.9 1.1],[mean_first(2) mean_last(2)],[sem_first(2) sem_last(2)],'k','Marker','o');
h.LineStyle = 'none';

for i = 1:num_subj
    plot([1.9 2.1],[NormThresh_first(i,3) NormThresh_last(i,3)], 'color',colormice(i,:),'Marker','o')
    hold on
    
end 

h=errorbar([1.9 2.1],[mean_first(3) mean_last(3)],[sem_first(3) sem_last(3)],'k','Marker','o');
h.LineStyle = 'none';

ylim([-0.02 0.12])
xlim([0.7 2.3])
hline([-15 15],'k:')
hline(0,'k:')
set(gca,'XTick',[1 2],'XTickLabel',{'15','-15'},'TickDir','out')
xlabel('Base Orien (Deg)')
ylabel('\Delta FA rate')