 Allmice = {'568','569','570','571'};

num_subj = size(Allmice,2);
Oriens = NaN(num_subj,2,7); % num of mice x led x orientations
OrienChoose = NaN(num_subj,1);
Threshall = NaN(num_subj,2);
Hit = NaN(num_subj,2,3); % num of mice x led x offs
FA = NaN(num_subj,2,3);
Thresh = NaN(num_subj,2,3);
RT_90 = NaN(num_subj,2,3);
RT_22 = NaN(num_subj,2,3);
FA_RT = NaN(num_subj,2,3); % only get the mean RT
Hit_90 = NaN(num_subj,2); % get the hit rate for 90 deg
Hit_22 = NaN(num_subj,2); % get the hit rate for 22.5 deg

for i_mice = 1:num_subj
    name = Allmice{i_mice};
    out = load(name);
    Oriens(i_mice,1,1:length(out.Output{1}.Infor.Orien)) = out.Output{1}.Infor.Orien;
    Oriens(i_mice,2,1:length(out.Output{2}.Infor.Orien)) = out.Output{2}.Infor.Orien;
    % chose the orientation in the second smallest one in the control
    % condition
    OrienChoose(i_mice,1) = out.Output{1}.Infor.Orien(2);
    for i_led = 1:2
        idx = [];
        idx = find(out.Output{i_led}.Infor.Orien == 90);
        Hit_90(i_mice,i_led)= out.Output{i_led}.target.all.HT_rate(idx); % get the hit rate for 90 deg
        idx = [];
        idx = intersect(find(out.Output{i_led}.Infor.Orien >= 20),find(out.Output{i_led}.Infor.Orien <= 29));
        Hit_22(i_mice,i_led) =  out.Output{i_led}.target.all.HT_rate(idx);; % get the hit rate for 22.5 deg
        
    end
    % get the hit rate and reaction time for that orientation
    for i_off = 1:3
        for i_led = 1:2
            idx = find(out.Output{i_led}.Infor.Orien == out.Output{1}.Infor.Orien(2));
            Hit(i_mice,i_led,i_off) = out.Output{i_led}.target.c_hit{i_off}(idx);
            idx = find(out.Output{i_led}.Infor.Orien == 90);
            RT_90(i_mice,i_led,i_off) = mean(out.Output{i_led}.target.RTonHit{i_off}{idx});
            idx = intersect(find(out.Output{i_led}.Infor.Orien >= 20),find(out.Output{i_led}.Infor.Orien <= 29));
            if ~isempty(idx)
            RT_22(i_mice,i_led,i_off) =mean(out.Output{i_led}.target.RTonHit{i_off}{idx});
            end
            FA_RT(i_mice,i_led,i_off) = mean(out.Output{i_led}.FA.FA_RT{i_off});
            
        end
    end
    
    Threshall(i_mice,:) = out.Threshall;
    Thresh(i_mice,:,:) = out.Thresh;
    FA(i_mice,1,:) = out.Output{1}.FA.FA;
    FA(i_mice,2,:) = out.Output{2}.FA.FA;
    
    
end 

%% plot the hit rate
figure
subplot(1,2,1)
ledcolor = {[0 0 0] [0.2 0.6 1]};
for i = 1:2
    temp = squeeze(Hit(:,i,:));
    h = errorbar([1 2 3], mean(temp,1),std(temp,[],1)./sqrt(num_subj),'Color',ledcolor{i});
    h.LineStyle = 'none';
    hold on
    scatter([1 2 3], mean(temp,1),'MarkerEdgeColor',ledcolor{i},'MarkerFaceColor',[1 1 1],'SizeData',40)
    
    
end 
xlim([0.5 3.5])
ylim([0 1])
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Hit rate')
% get the normalized value 
subplot(1,2,2)
for i = 1:2
    temp_old = squeeze(Hit(:,i,:));
    temp = temp_old./repmat(temp_old(:,1),1,3);
    h = errorbar([1 2 3], mean(temp,1),std(temp,[],1)./sqrt(num_subj),'Color',ledcolor{i});
    h.LineStyle = 'none';
    hold on
    scatter([1 2 3], mean(temp,1),'MarkerEdgeColor',ledcolor{i},'MarkerFaceColor',[1 1 1],'SizeData',40)
    
    
end 
xlim([0.5 3.5])
% ylim([0 1])
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Hit rate')
%% plot the threshold 
idx_all = logical([1 1 1 1]);
figure
subplot(2,2,1)
ledcolor = {[0 0 0] [0.2 0.6 1]};
mousecolor = {[0.5 0.5 0.5] [0.4 0.6 0.6]};
num = sum(idx_all);
for i = [2 1]
    temp = [];
    temp = squeeze(Thresh(idx_all,i,:));
    % plot the individual mouse
    for i_mice = 1:num
    scatter([1 2 3], temp(i_mice,:),'MarkerEdgeColor',mousecolor{i},'SizeData',10)
    hold on
    end 
    h = errorbar([1 2 3], mean(temp,1),std(temp,[],1)./sqrt(num),'Color',ledcolor{i});
    h.LineStyle = 'none';
    h.Marker = 'o';
   
    
end 
% statistic test, 2 way anova
% anovatest = zeros(8,3);
% anovatest(1:4,:) = squeeze(Thresh(idx_all,1,:));
% anovatest(5:8,:) = squeeze(Thresh(idx_all,2,:));
% 
% [p, tbl,stats] = anova2(anovatest,4);
% 
% %

text(1,40,['n = ', num2str(num),'mice'])
xlim([0.7 3.3])
ylim([0 50])
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Threshold (deg)')
% get the normalized value 
subplot(2,2,2)
for i = [2 1]
    temp_old = squeeze(Thresh(idx_all,i,:));
    temp = temp_old./repmat(temp_old(:,1),1,3);
    for i_mice = 1:num
    scatter([1 2 3], temp(i_mice,:),'MarkerEdgeColor',mousecolor{i},'SizeData',10)
    hold on
    end 
    
    h = errorbar([1 2 3], mean(temp,1),std(temp,[],1)./sqrt(num),'Color',ledcolor{i});
    h.LineStyle = 'none';
     h.Marker = 'o';
end 

% % % statistic test, 2 way anova
% anovatest = zeros(8,3);
% temp_old = squeeze(Thresh(idx_all,1,:));
% temp = temp_old./repmat(temp_old(:,1),1,3);
% anovatest(1:4,:) = temp;
% temp_old = squeeze(Thresh(idx_all,2,:));
% temp = temp_old./repmat(temp_old(:,1),1,3);
% anovatest(5:8,:) = temp;
% 
% [p, tbl,stats] = anova2(anovatest,4);

% [ht,pt]= ttest(anovatest(1:2,3),anovatest(3:4,3),'Tail','left');
% %
xlim([0.7 3.3])
ylim([0 1])
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Norm. Threshold')
subplot(2,2,3)
for i = [2 1]
    temp = [];
    temp = squeeze(FA(idx_all,i,:));
    for i_mice = 1:num
    scatter([1 2 3], temp(i_mice,:),'MarkerEdgeColor',mousecolor{i},'SizeData',10)
    hold on
    end 
    h = errorbar([1 2 3], mean(temp,1),std(temp,[],1)./sqrt(num),'Color',ledcolor{i});
    h.LineStyle = 'none';
    h.Marker = 'o';
    
end 
% % statistic test, 2 way anova
% anovatest = zeros(8,3);
% anovatest(1:4,:) = squeeze(FA(idx_all,1,:));
% anovatest(5:8,:) = squeeze(FA(idx_all,2,:));
% 
% [p, tbl,stats] = anova2(anovatest,4);
% 
% %

xlim([0.7 3.3])
ylim([0 0.3])
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('FA rate')
% get the normalized value 
subplot(2,2,4)
for i = [2 1]
    temp_old = squeeze(FA(idx_all,i,:));
    temp = temp_old./repmat(temp_old(:,1),1,3);
    
    for i_mice = 1:num
    scatter([1 2 3], temp(i_mice,:),'MarkerEdgeColor',mousecolor{i},'SizeData',10)
    hold on
    end 
    
    h = errorbar([1 2 3], mean(temp,1),std(temp,[],1)./sqrt(num),'Color',ledcolor{i});
    h.LineStyle = 'none';
   
     h.Marker = 'o';
    
end 

% % statistic test, 2 way anova
% anovatest = zeros(8,3);
% temp_old = squeeze(FA(idx_all,1,:));
% temp = temp_old./repmat(temp_old(:,1),1,3);
% anovatest(1:4,:) = temp;
% temp_old = squeeze(FA(idx_all,2,:));
% temp = temp_old./repmat(temp_old(:,1),1,3);
% anovatest(5:8,:) = temp;
% 
% [p, tbl,stats] = anova2(anovatest,4);
% 
% [ht,pt]= ttest(anovatest(1:2,3),anovatest(3:4,3),'Tail','Right');
% %

xlim([0.7 3.3])
ylim([0 4])
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('Norm. FA rate')



%% plot the RT difference between led and non led  

idx_all = logical([0 0 1 1 0]);
figure
subplot(1,3,1)
ledcolor = {[0 0 0] [0.2 0.6 1]};
mousecolor = {[0.5 0.5 0.5] [0.4 0.6 0.6]};
num = sum(idx_all);
for i = [2 1]
    temp = [];
    temp = squeeze(RT_90(idx_all,i,:))./1000;
    % plot the individual mouse
    for i_mice = 1:num
    scatter([1 2 3], temp(i_mice,:),'MarkerEdgeColor',mousecolor{i},'MarkerFaceColor',mousecolor{i},'SizeData',10)
    hold on
    end 
    h = errorbar([1 2 3], mean(temp,1),std(temp,[],1)./sqrt(num),'Color',ledcolor{i});
    h.LineStyle = 'none';
    
    hold on
    scatter([1 2 3], mean(temp,1),'MarkerEdgeColor',ledcolor{i},'MarkerFaceColor',[1 1 1],'SizeData',50)
    
    
end 
text(1,0.5,['n = ', num2str(num),'mice'])
xlim([0.7 3.3])
ylim([0 0.6])
hline([0.2 0.55], 'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('RT on Hits (90 deg)')

subplot(1,3,2)
ledcolor = {[0 0 0] [0.2 0.6 1]};
mousecolor = {[0.5 0.5 0.5] [0.4 0.6 0.6]};
num = sum(idx_all);
for i = [2 1]
    temp = [];
    temp = squeeze(RT_22(idx_all,i,:))./1000;
    % plot the individual mouse
    for i_mice = 1:num
    scatter([1 2 3], temp(i_mice,:),'MarkerEdgeColor',mousecolor{i},'MarkerFaceColor',mousecolor{i},'SizeData',10)
    hold on
    end 
    h = errorbar([1 2 3], mean(temp,1),std(temp,[],1)./sqrt(num),'Color',ledcolor{i});
    h.LineStyle = 'none';
    
    hold on
    scatter([1 2 3], mean(temp,1),'MarkerEdgeColor',ledcolor{i},'MarkerFaceColor',[1 1 1],'SizeData',50)
    
    
end 
text(1,0.5,['n = ', num2str(num),'mice'])
xlim([0.7 3.3])
ylim([0 0.6])
hline([0.2 0.55], 'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('RT on Hits (22.5 deg)')
subplot(1,3,3)
ledcolor = {[0 0 0] [0.2 0.6 1]};
mousecolor = {[0.5 0.5 0.5] [0.4 0.6 0.6]};
num = sum(idx_all);
for i = [2 1]
    temp = [];
    temp = squeeze(FA_RT(idx_all,i,:))./1000;
    % plot the individual mouse
    for i_mice = 1:num
    scatter([1 2 3], temp(i_mice,:),'MarkerEdgeColor',mousecolor{i},'MarkerFaceColor',mousecolor{i},'SizeData',10)
    hold on
    end 
    h = errorbar([1 2 3], mean(temp,1),std(temp,[],1)./sqrt(num),'Color',ledcolor{i});
    h.LineStyle = 'none';
    
    hold on
    scatter([1 2 3], mean(temp,1),'MarkerEdgeColor',ledcolor{i},'MarkerFaceColor',[1 1 1],'SizeData',50)
    
    
end 
text(1,0.5,['n = ', num2str(num),'mice'])
xlim([0.7 3.3])
ylim([0 0.6])
hline([0.2 0.55], 'k:')
set(gca,'XTick',[1 2 3],'XTickLabel',{'250','500','750'},'TickDir','out')
xlabel('ISI (ms)')
ylabel('RT on FA')