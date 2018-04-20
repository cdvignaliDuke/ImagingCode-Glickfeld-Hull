%% summary of FA rate of all mice
% need to wait until i558 get enought trials
% incoporate bootstrap into the bottom analysis code..
Allmice = {'i527','i544','i549','i557'};
fitS = {};
fitpre = {};
dobootstrap=0;
S_fit = {};
L_fit = {};
fit_bin = {};
Avg_pre_fit = {};
Avg_pre_thresh = zeros(3,3,size(Allmice,2));
Avg_pre_FA = zeros(3,3,size(Allmice,2));
% bootStats = {};
% bootStatspre = {};
% bootStats_sep = {};
FA_bin = zeros(size(Allmice,2),3);
skew = zeros(size(Allmice,2),4);
FA = zeros(size(Allmice,2),3);
S_FA = zeros(size(Allmice,2),3);
L_FA = zeros(size(Allmice,2),3);
Re = zeros(size(Allmice,2),3);
S_Re = zeros(size(Allmice,2),3);
L_Re = zeros(size(Allmice,2),3);
Thresh_Trl = [];
FA_Trl = [];
Length_Hit = zeros(size(Allmice,2),3);
Length_FA = zeros(size(Allmice,2),3);
MLength_Hit = zeros(size(Allmice,2),3);% matched length 
MLength_FA = zeros(size(Allmice,2),3); % matched length
MThresh =  zeros(size(Allmice,2),3); % matched threshold
MFA =  zeros(size(Allmice,2),3); % matched FA
Thresh = zeros(size(Allmice,2),3);
Threshpre = zeros(size(Allmice,2),3);
Thresh_S = zeros(size(Allmice,2),3);
Thresh_L = zeros(size(Allmice,2),3);
Thresh_bin = zeros(size(Allmice,2),3);
NBootstrapReps = 1000;
FA_pre = zeros(size(Allmice,2),3);
FA_pre_confi = zeros(3,2,size(Allmice,2));

FA_pre2 = zeros(size(Allmice,2),3);
FA_pre2_confi = zeros(3,2,size(Allmice,2));

FA_sep = zeros(3,3,size(Allmice,2));
FA_sep_confi = zeros(3,3,2,size(Allmice,2));

D_prime =  zeros(3,5,size(Allmice,2));
criterion = zeros(3,5,size(Allmice,2));

RTonHit =  zeros(3,5,size(Allmice,2));
RTonFA = zeros(3,size(Allmice,2));
Orien = zeros(size(Allmice,2),5);
FA_confi =  zeros(3,2,size(Allmice,2));


for i_mice = 1: size(Allmice,2)
    name = [Allmice{i_mice} '_' 'ISI-bin.mat'];
    load(name)
    % trial length
    Thresh_Trl(i_mice,:) = output.Trl.thresh;
    FA_Trl (i_mice,:) = output.Trl.FA;
    Length_Hit (i_mice,:) = output.Trl.length.Hit;
    Length_FA (i_mice,:) = output.Trl.length.FA;
    MLength_Hit (i_mice,:) = output.Trl.new.HM_length_mean;
    MLength_FA (i_mice,:) = output.Trl.new.FACR_lengthmean;
    MThresh(i_mice,:) =  output.Trl.new.thresh; % matched threshold
    MFA(i_mice,:) =  output.Trl.new.FA; % matched FA rate
    
    
    FA(i_mice,:)=output.FA.FA';
    FA_bin(i_mice,:) = output.FAbin.FA';
    FA_confi(1:3,1:2,i_mice) = output.FA.FA_confi;
    FA_pre(i_mice,:)=output.con.pre.FA';
    FA_pre_confi(1:3,1:2,i_mice) = output.con.pre.FA_confi;
    FA_pre2(i_mice,:)=output.con.pre2.FA';
    FA_pre2_confi(1:3,1:2,i_mice) = output.con.pre2.FA_confi;
    
    FA_bin_pre2(i_mice,:) = output.conbin.pre2.FA;
    
    FA_sep(:,:,i_mice) = output.con.FA;
    FA_sep_confi(:,:,:,i_mice) = output.con.FA_confi;
    
    S_FA(i_mice,:)=output.FA.S_FA';
    L_FA(i_mice,:)=output.FA.L_FA';
    Re(i_mice,:)=output.Re.Re';
    S_Re(i_mice,:)=output.Re.S_Re';
    L_Re(i_mice,:)=output.Re.L_Re';
    if length(output.Infor.Orien)==5
    Orien(i_mice,:)= output.Infor.Orien;
    
    D_prime(1:3,1:5,i_mice) = reshape(cell2mat(output.sdt.dprime),5,3)';
    criterion(1:3,1:5,i_mice) =  reshape(cell2mat(output.sdt.criterion),5,3)';
    temp = [];
    temp = cell2mat(output.target.RTonHit_mean);
    RTonHit(1:3,1:5,i_mice) = temp(1:2:5,:);
    RTonFA(1:3,i_mice) = cellfun(@mean, output.FA.FA_RT);
    else
    Orien(i_mice,1:4)= output.Infor.Orien;
    
    D_prime(1:3,1:4,i_mice) = reshape(cell2mat(output.sdt.dprime),4,3)';
    criterion(1:3,1:4,i_mice) =  reshape(cell2mat(output.sdt.criterion),4,3)';
    temp = [];
    temp = cell2mat(output.target.RTonHit_mean);
    RTonHit(1:3,1:4,i_mice) = temp(1:2:5,:);
    RTonFA(1:3,i_mice) = cellfun(@mean, output.FA.FA_RT);   
    end
    
    
    Ix_orien = find((Orien(i_mice,:)>15 & Orien(i_mice,:)<25)==1);
    %Ix_orien = 1;
    d_22Deg(i_mice,:) = squeeze(D_prime(:,Ix_orien,i_mice))';
    if length(output.Infor.Orien)==5
    d_90Deg(i_mice,:) = squeeze(D_prime(:,5,i_mice))';
    RT_22Deg(i_mice,:) = squeeze(RTonHit(:,Ix_orien,i_mice))';
    RT_90Deg(i_mice,:) = squeeze(RTonHit(:,5,i_mice))';
    
    c_22Deg(i_mice,:) = squeeze(criterion(:,Ix_orien,i_mice))';
    c_90Deg(i_mice,:) = squeeze(criterion(:,5,i_mice))';
    else
    d_90Deg(i_mice,:) = squeeze(D_prime(:,4,i_mice))';
    RT_22Deg(i_mice,:) = squeeze(RTonHit(:,Ix_orien,i_mice))';
    RT_90Deg(i_mice,:) = squeeze(RTonHit(:,4,i_mice))';
    
    c_22Deg(i_mice,:) = squeeze(criterion(:,Ix_orien,i_mice))';
    c_90Deg(i_mice,:) = squeeze(criterion(:,4,i_mice))';  
    end 
    
    
    
    % get the bootstrap working and include the error bar within mice
    for i_off = 1:length(output.Infor.Off)
        trialVec = (output.target.HT_num{i_off,1}+output.target.Miss_num{i_off,1})';
        trialVpre = (output.pre.HT_num{i_off,1}+output.pre.Miss_num{i_off,1})';
        fitS{i_mice,1}{i_off,1} = weibullFitLG(output.Infor.Orien, output.target.c_hit{i_off,1}',1, 1, {'nTrials', trialVec}); % use
        fitPre{i_mice,1}{i_off,1} = weibullFitLG(output.Infor.Orien, output.pre.c_hit{i_off,1}',1, 1, {'nTrials', trialVpre});
        Thresh(i_mice,i_off) = fitS{i_mice,1}{i_off,1}.thresh;
        Threshpre(i_mice,i_off) = fitPre{i_mice,1}{i_off,1}.thresh;
        
        trialVec = (output.target.S_HT_num{i_off,1}+output.target.S_Miss_num{i_off,1})';
        S_fit{i_mice,1}{i_off,1} = weibullFitLG(output.Infor.Orien, output.target.S_hit{i_off,1}',1, 1, {'nTrials', trialVec}); 
        
        trialVec = (output.target.L_HT_num{i_off,1}+output.target.L_Miss_num{i_off,1})';
        L_fit{i_mice,1}{i_off,1} = weibullFitLG(output.Infor.Orien, output.target.L_hit{i_off,1}',1, 1, {'nTrials', trialVec}); 
   
        Thresh_S(i_mice,i_off) =S_fit{i_mice,1}{i_off,1}.thresh;
        Thresh_L(i_mice,i_off) =L_fit{i_mice,1}{i_off,1}.thresh;
        
        trialVec = (output.targetbin.HT_num{i_off,1}+output.targetbin.Miss_num{i_off,1})';
        fit_bin{i_mice,1}{i_off,1} = weibullFitLG(output.Infor.Orien, output.targetbin.c_hit{i_off,1}',1, 1, {'nTrials', trialVec}); 
        
         
        Thresh_bin(i_mice,i_off) = fit_bin{i_mice,1}{i_off,1}.thresh;
        
        trialVec = (output.prebin.HT_num{i_off,1}+output.prebin.Miss_num{i_off,1})';
        fit_prebin{i_mice,1}{i_off,1} = weibullFitLG(output.Infor.Orien,output.prebin.c_hit{i_off,1}',1, 1, {'nTrials', trialVec}); 
        
         
        Thresh_prebin(i_mice,i_off) = fit_prebin{i_mice,1}{i_off,1}.thresh;
        
        
        for ii_off = 1:length(output.Infor.Off)
            
            Avg_pre_FA(i_off,ii_off,i_mice) = output.conbin.FA(i_off,ii_off); % first off is the pre pre off
            trialVec = ( squeeze(output.conbin.HT_num(i_off,ii_off,:))+ squeeze(output.conbin.Miss_num(i_off,ii_off,:)))';
            
            Avg_pre_fit{i_mice,1}{i_off,ii_off} = weibullFitLG(output.Infor.Orien, squeeze(output.conbin.c_hit(i_off,ii_off,:))',1, 1, {'nTrials', trialVec});
            
            
            Avg_pre_thresh(i_off, ii_off, i_mice) =  Avg_pre_fit{i_mice,1}{i_off,ii_off}.thresh;
        end
        hit_22Deg(i_mice,i_off) = output.target.c_hit{i_off,1}(Ix_orien);
        hit_22Deg_confi(i_mice,1:2,i_off) = output.target.confi{i_off,1}(Ix_orien,:);
        if length(output.Infor.Orien)==5
        
        hit_90Deg(i_mice,i_off) = output.target.c_hit{i_off,1}(5);
        else
        hit_90Deg(i_mice,i_off) = output.target.c_hit{i_off,1}(4);
        end
        % do bootstrap fitting
        
        skew(i_mice,i_off) = skewness(output.RT.B_RT{i_off}); 
        meanRT(i_mice,i_off) = mean(output.RT.B_RT{i_off});

        if dobootstrap==1
            [bootStats{i_mice,1}{i_off,1}]=BootstrapWeibullFit(trialVec, output.target.c_hit{i_off,1}',NBootstrapReps,output.Infor.Orien,1, 1);
            
            [bootStatspre{i_mice,1}{i_off,1}]=BootstrapWeibullFit(trialVpre,output.pre.c_hit{i_off,1}',NBootstrapReps,output.Infor.Orien,1, 1);
            
            for ii_off = 1:length(output.Infor.Off)
                trialVec_sep = (squeeze(output.con.HT_num(i_off,ii_off,:)) + squeeze(output.con.Miss_num(i_off,ii_off,:)))';
                [bootStats_sep{i_mice,1}{i_off,1}{ii_off,1}]=BootstrapWeibullFit(trialVec_sep,squeeze(output.con.c_hit(i_off,ii_off,:))',NBootstrapReps,output.Infor.Orien,1, 1);
                Thresh_sep(i_mice,i_off,ii_off) = bootStats_sep{i_mice,1}{i_off,1}{ii_off,1}.threshMean;
                
                
            end
        end
        
        % get the only condition when seperate by previous offs, and fixed in
        % distance of 750ms intervals, then
        % squeeze(output.con.c_hit(i_off,3,:))'
        
        
    end
     skew (i_mice,4) = skewness(output.RT.T_orien{1}(output.RT.T_orien{1}<=550  ));
     meanRT(i_mice,4) = mean(output.RT.T_orien{1}(output.RT.T_orien{1}<=550));
     temp = [];
     temp = output.RT.T_orien{1}(output.RT.T_orien{1}<=550 );
    
    
end
cbin_FA = output.Trl.cbin.FA;
cbin_Hit = output.Trl.cbin.Hit;


%% summary the skew plot 
figure

 scatter([250 500 750 1000], mean(skew,1),'MarkerEdgeColor',[0 0 0],'SizeData',80) 
 hold on
 for i = 1:size(skew,1)
     color_temp = [0.5 0.5 0.5];
     if i==2
         color_temp = [1 0 0];
     end
 scatter([250 500 750 1000],skew(i,:),'MarkerEdgeColor',color_temp,'SizeData',20)
 end
h= errorbar([250 500 750 1000], mean(skew,1),std(skew,[],1)./sqrt(size(skew,1)),'Color',[0 0 0]);
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'none';

  ylim([-1.5 1])
 xlim([200 1050])
  set(gca,'TickDir','out')
  set(gca, 'XTick',[250 500 750 1000])
  set(gca,'XTickLabel',{'250','500','750','22.5 deg'})
  ylabel('skewness')
  axis square

 
  
   [p,t,stats] = anova1(meanRT);
  
[c,m,h,nms] = multcompare(stats);

[h,p]=ttest(meanRT(:,3),meanRT(:,4));
%% plot the summayr RT on hit and RT on FA
Colors =flipud( gray(6));
figure
subplot(1,2,1)
a=[];
b=[];
for i_orien = 1:5
    scatter(squeeze(RTonHit(1,i_orien,:)),squeeze(RTonHit(3,i_orien,:)),'MarkerEdgeColor',Colors(i_orien+1,:),'SizeData',80)
    a=[a;squeeze(RTonHit(1,i_orien,:)) ];
    b= [b;squeeze(RTonHit(3,i_orien,:))];
    hold on
end 

plot([200 550], [200 550],'k:')
ylim([200 500])
xlim([200 500])
set(gca, 'TickDir','out')
[h,p]= ttest2(a,b);
xlabel('RT at Hit-250 (ms)')
ylabel('RT at Hit-750 (ms)')


subplot(1,2,2)

scatter(RTonFA(1,:),RTonFA(3,:),'MarkerEdgeColor',[0 0 0],'SizeData',80)
hold on
plot([200 550], [200 550],'k:')
ylim([200 500])
xlim([200 500])
set(gca, 'TickDir','out')
xlabel('RT at FA-250 (ms)')
ylabel('RT at FA-750 (ms)')
[h,p]= ttest2(RTonFA(1,:),RTonFA(3,:));
%% lot the short cycle vs long cycle on FA and normalized threshold 

subplot(1,2,1)
scatter(output.Infor.Off, mean(S_FA,1),'MarkerEdgeColor',[0 0 0],'SizeData',80)
hold on
h= errorbar(output.Infor.Off, mean(S_FA,1),std(S_FA,[],1)./sqrt(size(S_FA,1)),'Color',[0 0 0]);
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'none';

scatter(output.Infor.Off, mean(L_FA,1),'MarkerEdgeColor',[0.5 0.5 0.5],'SizeData',80)
hold on
h= errorbar(output.Infor.Off, mean(L_FA,1),std(L_FA,[],1)./sqrt(size(L_FA,1)),'Color',[0.5 0.5 0.5]);
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'none';


ylim([0 0.12])
xlim([200 800])
ylabel('FA rate ')
set(gca,'XTick',output.Infor.Off)
set(gca,'YTick',[0:0.03:0.12])
set(gca,'TickDir','out')
xlabel('ISI(ms)')
axis square

subplot(1,2,2)


Thresh_S_N = Thresh_S./repmat(Thresh_S(:,1),1,3); 
Thresh_L_N = Thresh_L./repmat(Thresh_L(:,1),1,3); 

scatter(output.Infor.Off, mean(Thresh_S_N,1),'MarkerEdgeColor',[0 0 0],'SizeData',80)
hold on
h= errorbar(output.Infor.Off, mean(Thresh_S_N,1),std(Thresh_S_N,[],1)./sqrt(size(Thresh_S_N,1)),'Color',[0 0 0]);
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'none';

scatter(output.Infor.Off, mean(Thresh_L_N,1),'MarkerEdgeColor',[0.5 0.5 0.5],'SizeData',80)
hold on
h= errorbar(output.Infor.Off, mean(Thresh_L_N,1),std(Thresh_L_N,[],1)./sqrt(size(Thresh_L_N,1)),'Color',[0.5 0.5 0.5]);
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'none';
ylim([0 1])
xlim([200 800])
ylabel('Tresh-black-short-gray long ')
set(gca,'XTick',output.Infor.Off)
set(gca,'YTick',[0:0.2:1])
set(gca,'TickDir','out')
xlabel('ISI(ms)')
axis square

% two way anova comparison
% [p,tbl,stats] = anova2([Thresh_S_N;Thresh_L_N],6)
% c= multcompare(stats,'Estimate','row');
%% plot the mean FA over bined mean FA
% figure
subplot(2,2,3)
scatter(output.Infor.Off, mean(FA,1),'MarkerEdgeColor',[0 0 0],'SizeData',80)
hold on
h= errorbar(output.Infor.Off, mean(FA,1),std(FA,[],1)./sqrt(size(FA,1)),'Color',[0 0 0]);
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'none';

scatter(output.Infor.Off, mean(FA_bin,1),'MarkerEdgeColor',[0.5 0.5 0.5],'SizeData',80)
hold on
h= errorbar(output.Infor.Off, mean(FA_bin,1),std(FA_bin,[],1)./sqrt(size(FA_bin,1)),'Color',[0.5 0.5 0.5]);
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'none';

ylim([0 0.12])
xlim([200 800])
ylabel('FA rate ')
set(gca,'XTick',output.Infor.Off)
set(gca,'YTick',[0:0.03:0.12])
set(gca,'TickDir','out')
xlabel('ISI(ms)')

subplot(2,2,4)

Thresh_N = Thresh./repmat(Thresh(:,1),1,3); 
Thresh_bin_N = Thresh_bin./repmat(Thresh(:,1),1,3); 

scatter(output.Infor.Off, mean(Thresh_N,1),'MarkerEdgeColor',[0 0 0],'SizeData',80)
hold on
h= errorbar(output.Infor.Off, mean(Thresh_N,1),std(Thresh_N,[],1)./sqrt(size(Thresh_N,1)),'Color',[0 0 0]);
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'none';

scatter(output.Infor.Off, mean(Thresh_bin_N,1),'MarkerEdgeColor',[0.5 0.5 0.5],'SizeData',80)
hold on
h= errorbar(output.Infor.Off, mean(Thresh_bin_N,1),std(Thresh_bin_N,[],1)./sqrt(size(Thresh_bin_N,1)),'Color',[0.5 0.5 0.5]);
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'none';
ylim([0 1.4])
xlim([200 800])
ylabel('Tresh ')
set(gca,'XTick',output.Infor.Off)
set(gca,'YTick',[0:0.2:1.4])
set(gca,'TickDir','out')
xlabel('ISI(ms)')

supertitle('black: n-off; gray: fixed offs')




%% plot summary FA rate
figure
% Colors = lines(16);

for i_mice = 1:size(Allmice,2)
    scatter(output.Infor.Off,FA(i_mice,:),'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'SizeData',20)
    hold on
    %scatter(output.Infor.Off,FA_pre2(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:))
%     for i_off = 1:length(output.Infor.Off)
%         line(repmat(output.Infor.Off(i_off),1,2),FA_confi(i_off,:,i_mice),'color',Colors(i_mice,:))
%     end
%     text(220,0.15-i_mice./100,Allmice{i_mice},'color',Colors(i_mice,:),'FontSize',14)
    ylim([0 0.12])
    xlim([200 800])
end
hold on
h=errorbar(output.Infor.Off,mean(FA,1),std(FA,[],1)./sqrt(size(Allmice,2)),'k');
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'none';
hold on
scatter(output.Infor.Off,mean(FA,1),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'SizeData',60)
% h=errorbar(output.Infor.Off,mean(FA_pre2,1),std(FA_pre2,[],1)./sqrt(size(Allmice,2)),'color',[0.5 0.5 0.5]);
% h.LineStyle = 'none';
% h.LineWidth = 1.3;
% h.Marker = 'o';


ylabel('FA rate ')
set(gca,'XTick',output.Infor.Off)
set(gca,'YTick',[0:0.03:0.12],'TickDir','out')
xlabel('ISI(ms)')
axis square


 [p1,tbl,stats] = anova1(FA);
 b= multcompare(stats);

%% plot FA rate separate by n offs, all n-1 trials
figure
Colors = lines(16);
subplot(2,2,1)
for i_mice = 1:size(Allmice,2)
    scatter(output.Infor.Off,FA_pre(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:),'MarkerFaceColor',Colors(i_mice,:))
    hold on
    %scatter(output.Infor.Off,FA_pre2(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:))
    for i_off = 1:length(output.Infor.Off)
        line(repmat(output.Infor.Off(i_off),1,2),FA_pre_confi(i_off,:,i_mice),'color',Colors(i_mice,:))
    end
    text(220,0.15-i_mice*2./100,Allmice{i_mice},'color',Colors(i_mice,:),'FontSize',10)
    ylim([0 0.2])
    xlim([200 800])
end
hold on
h=errorbar(output.Infor.Off,mean(FA_pre,1),std(FA_pre,[],1)./sqrt(size(Allmice,2)),'k');
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'o';
ylabel('FA rate ')
set(gca,'XTick',output.Infor.Off)
set(gca,'YTick',[0:0.05:0.2])
xlabel('ISI(ms)- n offs')
title('FA rate - All n-1 trials')

for i_off = 1:3  % different n-1 offs
    subplot(2,2,i_off+1)
    for i_mice = 1:size(Allmice,2)
        scatter(output.Infor.Off,squeeze(FA_sep(i_off,:,i_mice)),'MarkerEdgeColor',Colors(i_mice,:),'MarkerFaceColor',Colors(i_mice,:))
        hold on
        %scatter(output.Infor.Off,FA_pre2(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:))
        for ii_off = 1:length(output.Infor.Off)
            line(repmat(output.Infor.Off(ii_off),1,2),squeeze(FA_sep_confi(i_off,ii_off,:,i_mice)),'color',Colors(i_mice,:))
        end
        text(220,0.15-i_mice*2./100,Allmice{i_mice},'color',Colors(i_mice,:),'FontSize',10)
        ylim([0 0.2])
        xlim([200 800])
    end
    hold on
    h=errorbar(output.Infor.Off,mean(squeeze(FA_sep(i_off,:,:)),2),std(squeeze(FA_sep(i_off,:,:)),[],2)./sqrt(size(Allmice,2)),'k');
    h.LineStyle = 'none';
    h.LineWidth = 1.3;
    h.Marker = 'o';
    ylabel('FA rate ')
    set(gca,'XTick',output.Infor.Off)
    set(gca,'YTick',[0:0.05:0.2])
    xlabel('ISI(ms)- n offs')
    title(['FA rate-' num2str(output.Infor.Off(i_off))  'n-1 trials'])
end

%% plot FA rate separate by n-1 offs, all n trials
figure
Colors = lines(16);
subplot(2,2,1)
for i_mice = 1:size(Allmice,2)
    scatter(output.Infor.Off,FA_pre2(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:),'MarkerFaceColor',Colors(i_mice,:))
    hold on
    %scatter(output.Infor.Off,FA_pre2(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:))
    for i_off = 1:length(output.Infor.Off)
        line(repmat(output.Infor.Off(i_off),1,2),FA_pre2_confi(i_off,:,i_mice),'color',Colors(i_mice,:))
    end
    text(220,0.2-i_mice*2./100,Allmice{i_mice},'color',Colors(i_mice,:),'FontSize',10)
    ylim([0 0.2])
    xlim([200 800])
end
hold on
h=errorbar(output.Infor.Off,mean(FA_pre2,1),std(FA_pre2,[],1)./sqrt(size(Allmice,2)),'k');
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'o';
ylabel('FA rate ')
set(gca,'XTick',output.Infor.Off)
set(gca,'YTick',[0:0.05:0.2])
xlabel('ISI(ms)- n-1 offs')
title('FA rate - All n trials')

for i_off = 1:3  % different n-1 offs
    subplot(2,2,i_off+1)
    for i_mice = 1:size(Allmice,2)
        scatter(output.Infor.Off,squeeze(FA_sep(:,i_off,i_mice)),'MarkerEdgeColor',Colors(i_mice,:),'MarkerFaceColor',Colors(i_mice,:))
        hold on
        %scatter(output.Infor.Off,FA_pre2(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:))
        for ii_off = 1:length(output.Infor.Off)
            line(repmat(output.Infor.Off(ii_off),1,2),squeeze(FA_sep_confi(ii_off,i_off,:,i_mice)),'color',Colors(i_mice,:))
        end
        %text(220,0.15-i_mice*2./100,Allmice{i_mice},'color',Colors(i_mice,:),'FontSize',10)
        ylim([0 0.2])
        xlim([200 800])
    end
    hold on
    h=errorbar(output.Infor.Off,mean(squeeze(FA_sep(:,i_off,:)),2),std(squeeze(FA_sep(:,i_off,:)),[],2)./sqrt(size(Allmice,2)),'k');
    h.LineStyle = 'none';
    h.LineWidth = 1.3;
    h.Marker = 'o';
    ylabel('FA rate ')
    set(gca,'XTick',output.Infor.Off)
    set(gca,'YTick',[0:0.05:0.2])
    xlabel('ISI(ms)- n-1 offs')
    title(['FA rate-' num2str(output.Infor.Off(i_off))  'n trials'])
end

%% plot FA rate and threshold when fixed collapsed all N offs
figure
subplot(1,2,1)
T_preratio = Threshpre./repmat(Threshpre(:,1),1,3);
    for i_mice = 1:size(Allmice,2)
        scatter(output.Infor.Off,T_preratio(i_mice,:),'MarkerEdgeColor',[0.5 0.5 0.5],'SizeData',20);
       hold on
        ylim([0 1.2])
        xlim([200 800])
    end
    hold on
    h=errorbar(output.Infor.Off,mean(T_preratio,1),std(T_preratio,[],1)./sqrt(size(Allmice,2)),'k');
    h.LineStyle = 'none';
    h.LineWidth = 1.3;
    h.Marker = 'none';
    hold on
    scatter(output.Infor.Off,mean(T_preratio,1),'MarkerEdgeColor',[0 0 0],'SizeData',60)
    hline(1,'k:')
    ylabel('Orien Thresh Ratio')
    set(gca,'XTick',output.Infor.Off)
    set(gca,'YTick',[0:0.4:1.2])
    set(gca,'TickDir','out')
    xlabel('ISI(ms)-n-1 offs')
    title(['Threshhold- all N trials'])

%  [p1,tbl,stats] = anova1(T_preratio);
%  b= multcompare(stats);


subplot(1,2,2)
for i_mice = 1:size(Allmice,2)
        scatter(output.Infor.Off,FA_pre2(i_mice,:),'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Sizedata',20)
        hold on
       
        ylim([0 0.12])
        xlim([200 800])
    end
    hold on
    h=errorbar(output.Infor.Off,mean(FA_pre2,1),std(FA_pre2,[],1)./sqrt(size(Allmice,2)),'k');
    h.LineStyle = 'none';
    h.LineWidth = 1.3;
    h.Marker = 'none';
    
    scatter(output.Infor.Off,mean(FA_pre2,1),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'Sizedata',60)
    
    ylabel('FA rate ')
    set(gca,'XTick',output.Infor.Off)
    set(gca,'YTick',[0:0.03:0.12])
    set(gca,'TickDir','out')
    xlabel('ISI(ms)- n-1 offs')
    title(['FA rate-all N trials'])
%% plot FA rate when fixed N-1: 750ms 
figure
i_off = 3;
for i_mice = 1:size(Allmice,2)
        scatter(output.Infor.Off,squeeze(FA_sep(:,i_off,i_mice)),'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Sizedata',20)
        hold on
       
        ylim([0 0.12])
        xlim([200 800])
    end
    hold on
    h=errorbar(output.Infor.Off,mean(squeeze(FA_sep(:,i_off,:)),2),std(squeeze(FA_sep(:,i_off,:)),[],2)./sqrt(size(Allmice,2)),'k');
    h.LineStyle = 'none';
    h.LineWidth = 1.3;
    h.Marker = 'none';
    
    scatter(output.Infor.Off,mean(squeeze(FA_sep(:,i_off,:)),2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'Sizedata',60)
    
    ylabel('FA rate ')
    set(gca,'XTick',output.Infor.Off)
    set(gca,'YTick',[0:0.06:0.12])
    set(gca,'TickDir','out')
    xlabel('ISI(ms)- n-1 offs')
    title(['FA rate-' num2str(output.Infor.Off(i_off))  'n trials'])
%% plot normalized FA rate
figure
Colors = lines(16);

for i_mice = 1:size(Allmice,2)
    scatter(output.Infor.Off,FA(i_mice,:)./FA(i_mice,1),'MarkerEdgeColor',Colors(i_mice,:),'MarkerFaceColor',Colors(i_mice,:))
    hold on
    %  scatter(output.Infor.Off,FA_pre2(i_mice,:)./FA_pre2(i_mice,1),'MarkerEdgeColor',Colors(i_mice,:))
    
    text(220,15-i_mice,Allmice{i_mice},'color',Colors(i_mice,:),'FontSize',14)
    ylim([0 15])
    xlim([200 800])
end
hold on
h=errorbar(output.Infor.Off,mean(FA./repmat(FA(:,1),1,3),1),std(FA./repmat(FA(:,1),1,3),[],1)./sqrt(size(Allmice,2)),'k');
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'o';
% h=errorbar(output.Infor.Off,mean(FA_pre2./repmat(FA_pre2(:,1),1,3),1),std(FA_pre2./repmat(FA_pre2(:,1),1,3),[],1)./sqrt(size(Allmice,2)),'color',[0.5 0.5 0.5]);
% h.LineStyle = 'none';
% h.LineWidth = 1.3;
% h.Marker = 'o';


ylabel('FA rate normalized to 250ms')
set(gca,'XTick',output.Infor.Off)
ylim([0 5])
xlabel('ISI(ms)')

%% plot the change of d prime, RTonHit and criterion for 90 degree change and near 22.5 degree



% plot the change of d prime over 250ms condition
% figure
% h=plot(output.Infor.Off,(d_22Deg-repmat(d_22Deg(:,1),1,3)));
% for  i_mice = 1: size(Allmice,2)
%     h(i_mice).Color = Colors(i_mice,:)
%  end
% hold on
% h=plot(output.Infor.Off,(d_90Deg-repmat(d_90Deg(:,1),1,3)));
% for  i_mice = 1: size(Allmice,2)
%     h(i_mice).Color = Colors(i_mice,:)
%  end

% plot the raw d prime value for 90 degree and 22.5 degree
figure
subplot(1,2,1)
h=plot(output.Infor.Off,d_22Deg);
for  i_mice = 1: size(Allmice,2)
    h(i_mice).Color =[0.5 0.5 0.5];% Colors(i_mice,:);
    h(i_mice).LineStyle ='--';
end
hold on
errorbar(output.Infor.Off,mean(d_22Deg,1),std(d_22Deg,[],1)./sqrt(size(Allmice,2)),'k--')
scatter(output.Infor.Off,mean(d_22Deg,1),'k')
h=plot(output.Infor.Off,d_90Deg);
for  i_mice = 1: size(Allmice,2)
    h(i_mice).Color = [0.5 0.5 0.5];%Colors(i_mice,:);
end
errorbar(output.Infor.Off,mean(d_90Deg,1),std(d_90Deg,[],1)./sqrt(size(Allmice,2)),'k')
scatter(output.Infor.Off,mean(d_90Deg,1),'k','fiiled')
ylim([0 5])
xlim([200 800])
ylabel('d prime')
set(gca,'XTick',output.Infor.Off)
xlabel('ISI(ms)')
title('fiiled:90deg; open:~22.5 deg')


subplot(1,2,2)
h=plot(output.Infor.Off,c_22Deg);
for  i_mice = 1: size(Allmice,2)
    h(i_mice).Color =[0.5 0.5 0.5];% Colors(i_mice,:)
    h(i_mice).LineStyle ='--';
end
hold on
errorbar(output.Infor.Off,mean(c_22Deg,1),std(c_22Deg,[],1)./sqrt(size(Allmice,2)),'k--')
scatter(output.Infor.Off,mean(c_22Deg,1),'k')
h=plot(output.Infor.Off,c_90Deg);
for  i_mice = 1: size(Allmice,2)
    h(i_mice).Color = [0.5 0.5 0.5];%Colors(i_mice,:)
end
errorbar(output.Infor.Off,mean(c_90Deg,1),std(c_90Deg,[],1)./sqrt(size(Allmice,2)),'k')
scatter(output.Infor.Off,mean(c_90Deg,1),'k','fiiled')
ylim([-1 2])
xlim([200 800])
ylabel('criterion')
set(gca,'XTick',output.Infor.Off)
xlabel('ISI(ms)')
title('fiiled:90deg; open:~22.5 deg')

% subplot(1,3,1)
% h=plot(output.Infor.Off,RT_22Deg);
% for  i_mice = 1: size(Allmice,2)
%     h(i_mice).Color = Colors(i_mice,:)
%     h(i_mice).LineStyle ='--';
% end
% hold on
% errorbar(output.Infor.Off,mean(RT_22Deg,1),std(RT_22Deg,[],1)./sqrt(size(Allmice,2)),'k--')
% scatter(output.Infor.Off,mean(RT_22Deg,1),'k')
% h=plot(output.Infor.Off,RT_90Deg);
% for  i_mice = 1: size(Allmice,2)
%     h(i_mice).Color = Colors(i_mice,:)
% end
% errorbar(output.Infor.Off,mean(RT_90Deg,1),std(RT_90Deg,[],1)./sqrt(size(Allmice,2)),'k')
% scatter(output.Infor.Off,mean(RT_90Deg,1),'k','fiiled')
% ylim([100 400])
% xlim([200 800])
% ylabel('RT on Hits')
% set(gca,'XTick',output.Infor.Off)
% xlabel('ISI(ms)')
% title('fiiled:90deg; open:~22.5 deg')
%% plot what proportion of hit rate change are due to either criterion change or sensitivity change
% calculete the propotion of Hit rate change are solely due to criterion,
% only focus on the biggest change in hit rate condition, 22.5 deg change



D_hit =hit_22Deg- repmat(hit_22Deg(:,1),1,3);
D_hit = D_hit(:,2:3);
hit_250 = hit_22Deg(:,1);
d_250 = d_22Deg(:,1);
c_250 = c_22Deg(:,1);

P_c = (normcdf(repmat(d_250./2,1,2) - c_22Deg(:,2:3)) - repmat(hit_250,1,2))./D_hit ; 
P_d = (normcdf(d_22Deg(:,2:3)./2- repmat(c_250,1,2)) - repmat(hit_250,1,2))./D_hit ; 

for i_mice = 1: size(Allmice,2)
    for i_off = 1:2 % only for 500ms off and 750ms off condition
        % for the fix c condition, varying d'
        temp = sort([d_22Deg(i_mice,1) d_22Deg(i_mice,i_off+1)]);
        x_temp =  fminbnd(@(x) (normcdf(x./2 - c_22Deg(i_mice,i_off+1)) - normcdf(x./2- c_22Deg(i_mice,1))), temp(1),temp(2));
        x_temp_max = fminbnd(@(x) -(normcdf(x./2 - c_22Deg(i_mice,i_off+1)) - normcdf(x./2- c_22Deg(i_mice,1))), temp(1),temp(2));
        P_c_min(i_mice,i_off) =(normcdf(x_temp./2 - c_22Deg(i_mice,i_off+1)) - normcdf(x_temp./2- c_22Deg(i_mice,1)))./D_hit(i_mice,i_off);
        P_c_max(i_mice,i_off) =(normcdf(x_temp_max./2 - c_22Deg(i_mice,i_off+1)) - normcdf(x_temp_max./2- c_22Deg(i_mice,1)))./D_hit(i_mice,i_off);
        
        % for the fix d condition, varying c
        temp = sort([c_22Deg(i_mice,1) c_22Deg(i_mice,i_off+1)]);
        x_temp =  fminbnd(@(x) (normcdf(d_22Deg(i_mice,i_off+1)./2 - x) - normcdf(d_22Deg(i_mice,1)./2- x)), temp(1),temp(2));
        x_temp_max = fminbnd(@(x) -(normcdf(d_22Deg(i_mice,i_off+1)./2 - x) - normcdf(d_22Deg(i_mice,1)./2- x)), temp(1),temp(2));
        P_d_min(i_mice,i_off) =(normcdf(d_22Deg(i_mice,i_off+1)./2 - x_temp) - normcdf(d_22Deg(i_mice,1)./2- x_temp))./D_hit(i_mice,i_off);
        P_d_max(i_mice,i_off) =(normcdf(d_22Deg(i_mice,i_off+1)./2 - x_temp_max) - normcdf(d_22Deg(i_mice,1)./2- x_temp_max))./D_hit(i_mice,i_off);
        
       
    end
    
end

figure
% for explained by c only 
 subplot(1,2,1)
h1=bar([400, 650], mean(P_c,1));
h1.BarWidth = 0.4;
% h1.FaceColor = [1,0.8,0.4];
% h1.EdgeColor = [1,0.8,0.4];
 h1.FaceColor = [0.5,0.5,0.5];
 h1.EdgeColor = [0.5,0.5,0.5];

hold on 
h2=errorbar ([400, 650],mean(P_c,1),std(P_c,[],1)./sqrt(size(P_c,1)));
h2.LineStyle = 'none';
h2.Color = [0.5,0.5,0.5];%[1,0.8,0.4];
h2.LineWidth = 1;

% for explained by the d' only 
h3=bar([500, 750], mean(P_d,1));
h3.BarWidth = 0.4;
% h3.FaceColor = [0.6,0.6,0.8];
% h3.EdgeColor = [0.6,0.6,0.8];
h3.FaceColor = [0,0,0];
h3.EdgeColor = [0,0,0];
hold on 
h4=errorbar ([500, 750],mean(P_d,1),std(P_d,[],1)./sqrt(size(P_d,1)));
h4.LineStyle = 'none';
h4.Color = [0,0,0];%[0.6,0.6,0.8];
h4.LineWidth =1;
legend([h1,h3],['\Delta' 'c'],['\Delta' 'd'''],'Location','Northwest')
legend boxoff
xlim([300 850])
ylim([0 1.2])
set(gca,'XTick',[450 700],'XTickLabel',{'500ms vs. 250ms' '750ms vs. 250ms'})
set(gca,'TickDir','out')
ylabel('Proportion')
title(['Proportion of ' '\Delta' 'Hit rate due to ' '\Delta' 'c or ' '\Delta' 'd'' alone'])

% plot the maximum and minimum range of change of delta d and delta c
% method:http://dx.doi.org/10.1016/j.neuron.2015.05.007


subplot(1,2,2)
% for max and min proportion of fix c, vary on d'
h1=bar([400, 650], mean(P_c_min,1));
h1.BarWidth = 0.15;
h1.FaceColor = [1,0.8,0.4];
h1.EdgeColor = [1,0.8,0.4];
hold on 
h2=errorbar ([400, 650],mean(P_c_min,1),std(P_c_min,[],1)./sqrt(size(P_c_min,1)));
h2.LineStyle = 'none';
h2.Color = [1,0.8,0.4];
h2.LineWidth = 1;
hold on
h3=bar([440, 690], mean(P_c_max,1));
h3.BarWidth = 0.15;
h3.FaceColor = [1,0.8,0.4].*0.8;
h3.EdgeColor = [1,0.8,0.4].*0.8;
hold on 
h4=errorbar ([440, 690],mean(P_c_max,1),std(P_c_max,[],1)./sqrt(size(P_c_max,1)));
h4.LineStyle = 'none';
h4.Color = [1,0.8,0.4].*0.8;
h4.LineWidth = 1;
hold on
% for fix in d' vary on c 
h5=bar([500, 750], mean(P_d_min,1));
h5.BarWidth = 0.15;
h5.FaceColor = [0.6,0.6,0.8];
h5.EdgeColor = [0.6,0.6,0.8];
hold on 
h6=errorbar ([500, 750],mean(P_d_min,1),std(P_d_min,[],1)./sqrt(size(P_d_min,1)));
h6.LineStyle = 'none';
h6.Color = [0.6,0.6,0.8];
h6.LineWidth =1;
hold on
h7=bar([540, 790], mean(P_d_max,1));
h7.BarWidth = 0.15;
h7.FaceColor = [0.6,0.6,0.8].*0.8;
h7.EdgeColor = [0.6,0.6,0.8].*0.8;
hold on 
h8=errorbar ([540, 790],mean(P_d_max,1),std(P_d_max,[],1)./sqrt(size(P_d_max,1)));
h8.LineStyle = 'none';
h8.Color = [0.6,0.6,0.8].*0.8;
h8.LineWidth =1;

legend([h1,h3,h5,h7],['\Delta' 'c' ' min proportion'],['\Delta' 'c' ' max proportion'],['\Delta' 'd''' ' min proportion'],['\Delta' 'd''' ' max proportion'],'Location','Northwest')
legend boxoff
xlim([300 850])
set(gca,'XTick',[450 700],'XTickLabel',{'500ms vs. 250ms' '750ms vs. 250ms'})
ylabel('Proportion')
title(['Theoretically range of portion of ' '\Delta' 'Hit rate'])
    
%% instead plot point graph
figure

for i = 1:size(P_c,1)
scatter([400, 650],P_c(i,:),'MarkerFaceColor',[0.5,0.5,0.5],'MarkerEdgeColor',[0.5,0.5,0.5],'SizeData',20);
hold on
end
hold on
h1=scatter([400, 650], mean(P_c,1),'MarkerEdgeColor',[0.5,0.5,0.5],'SizeData',60);


hold on 
h2=errorbar ([400, 650],mean(P_c,1),std(P_c,[],1)./sqrt(size(P_c,1)));
h2.LineStyle = 'none';
h2.Color = [0.5,0.5,0.5];%[1,0.8,0.4];
h2.LineWidth = 1;

% for explained by the d' only 
for i = 1:size(P_d,1)
scatter([500, 750],P_d(i,:),'MarkerFaceColor',[0,0,0],'MarkerEdgeColor',[0,0,0],'SizeData',20);
end

h3=scatter([500, 750], mean(P_d,1),'MarkerEdgeColor',[0,0,0],'SizeData',60);



h4=errorbar ([500, 750],mean(P_d,1),std(P_d,[],1)./sqrt(size(P_d,1)));
h4.LineStyle = 'none';
h4.Color = [0,0,0];%[0.6,0.6,0.8];
h4.LineWidth =1;
legend([h1,h3],['\Delta' 'c'],['\Delta' 'd'''],'Location','Northwest')
legend boxoff
xlim([300 850])
ylim([-0.08 1.2])
set(gca,'XTick',[450 700],'XTickLabel',{'500ms vs. 250ms' '750ms vs. 250ms'})
set(gca,'TickDir','out','YTick',[0:0.3:1.2])
ylabel('Proportion')
title(['Proportion of ' '\Delta' 'Hit rate due to ' '\Delta' 'c or ' '\Delta' 'd'' alone'])
% use two way anova
[p,tbl,stats] = anova2([P_c;P_d],6)
c= multcompare(stats,'Estimate','row');

%% plot the ROC curve and the hit and FA rate fall into the curve

Colors = lines(16); % for the mice
colorchoice = {[0 0 0] [0.5 0.5 0.5] [1 0.8 0.4] [1 0 0]}; %  colors between different offs

figure

a=linspace(0,1,100);
for c=[0 0.5 1 1.5]
    b = normcdf(-2*c-norminv(a));
    plot(a,b,'k:')
    hold on
end
for d=[1 1.5 2]
    b = normcdf(d+norminv(a));
    plot(a,b,'k:')
    hold on
end
% now plot the animals's behavior on different sessions 
for i_mice = 1:size(Allmice,2)
    for i_off = 1:3
    scatter(FA(i_mice,i_off),hit_22Deg(i_mice,i_off),'MarkerEdgeColor',colorchoice{i_off},'MarkerFaceColor',colorchoice{i_off})
    hold on
    plot(FA(i_mice,i_off)*[1 1],squeeze(hit_22Deg_confi(i_mice,1:2,i_off)),'Color',colorchoice{i_off})
    plot(squeeze(FA_confi(i_off,1:2,i_mice)),hit_22Deg(i_mice,i_off)*[1 1],'Color',colorchoice{i_off})
    end
    hold on
    plot(FA(i_mice,:),hit_22Deg(i_mice,:),'Color',Colors(i_mice,:))
end

xlim([0 0.12])
ylim([0 1])
ylabel('Hit rate')
xlabel('FA rate')
title('22.5 Deg-c: [0 0.5 1 1.5] d:[1 1.5 2]')
%% plot the FA for long vs short cycle
figure
h1=errorbar(output.Infor.Off,mean(S_FA,1).*100,std(S_FA.*100,[],1)./sqrt(size(Allmice,2)),'color',[1 0.4 0.2]);
h1.LineStyle = 'none';
h1.LineWidth = 0.5;
h1.Marker = 'o';

hold on
h1=errorbar(output.Infor.Off,mean(L_FA,1).*100,std(L_FA.*100,[],1)./sqrt(size(Allmice,2)),'color',[0 0.6 0.4]);
h1.LineStyle = 'none';
h1.LineWidth = 0.5;
h1.Marker = 'o';

hold on
h=errorbar(output.Infor.Off,mean(FA,1).*100,std(FA.*100,[],1)./sqrt(size(Allmice,2)),'k');
h.LineStyle = 'none';
h.LineWidth = 1.1;
h.Marker = 'o';

legend('Cycle:3-5','Cycle:7-9','ALL','Location','northwest')
ylim([0 15])
xlim([200 800])
ylabel('FA rate (%)')
set(gca,'XTick',output.Infor.Off)
set(gca,'YTick',[0:3:15])
xlabel('ISI(ms)')



%% plot the threshold

figure
for i_mice = 1:size(Allmice,2)
    scatter(output.Infor.Off,Thresh(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:),'MarkerFaceColor',Colors(i_mice,:))
    hold on
    %  scatter(output.Infor.Off,Threshpre(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:))
    text(220,15-i_mice*2,Allmice{i_mice},'color',Colors(i_mice,:),'FontSize',14)
    for i_off = 1:3
        plot(output.Infor.Off(i_off)*[1 1], bootStats{i_mice,1}{i_off,1}.ci95,'Color',Colors(i_mice,:) );
        hold on
    end
    ylim([0 40])
    xlim([200 800])
end
hold on
h=errorbar(output.Infor.Off,mean(Thresh,1),std(Thresh,[],1)./sqrt(size(Allmice,2)),'k');
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'o';
% hold on
% h=errorbar(output.Infor.Off,mean(Threshpre,1),std(Threshpre,[],1)./sqrt(size(Allmice,2)),'r');
% h.LineStyle = 'none';
% h.LineWidth = 1.3;
% h.Marker = 'o';

ylabel('Orien Thresh (deg)')
set(gca,'XTick',output.Infor.Off)
% set(gca,'YTick',[0:3:15])
xlabel('ISI(ms)')

%  [p1,tbl,stats] = anova1(Thresh);
%  b= multcompare(stats);

%% plot thresh ratio relative to 250 ms off
T_ratio = Thresh./repmat(Thresh(:,1),1,3);
T_preratio = Threshpre./repmat(Threshpre(:,1),1,3);

figure
for i_mice = 1:size(Allmice,2)
    scatter(output.Infor.Off,T_ratio(i_mice,:),'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'SizeData',20);
   
    hold on
    %  h= scatter(output.Infor.Off,T_preratio(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:));
    %  h.SizeData = 55;
%     text(220,i_mice*0.1,Allmice{i_mice},'color',Colors(i_mice,:),'FontSize',14)
    ylim([0 1])
    xlim([200 800])
end
hold on
h=errorbar(output.Infor.Off,mean(T_ratio,1),std(T_ratio,[],1)./sqrt(size(Allmice,2)),'k');
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'none';
hold on
scatter(output.Infor.Off,mean(T_ratio,1),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'SizeData',60)
% h=errorbar(output.Infor.Off,mean(T_preratio,1),std(T_preratio,[],1)./sqrt(size(Allmice,2)),'color',[0.5 0.5 0.5]);
% h.LineStyle = 'none';
% h.LineWidth = 1.3;
% h.Marker = 'o';

ylabel('Orien Thresh Ratio (over ISI 250 ms)')
set(gca,'XTick',output.Infor.Off)
% set(gca,'YTick',[0:3:15])
xlabel('ISI(ms)')

   [p1,tbl,stats] = anova1(T_ratio);
   b= multcompare(stats);

%% plot the threshold difference between different offs
figure
h={};
colorchoice = {[0 0 0] [0.5 0.5 0.5] [1 0.8 0.4] [1 0 0]};
for  i_mice = 1:size(Allmice,2)
    subplot(2,2,i_mice)
    for i_off =1:3 %[1,3]
        
        maxI = max(output.Infor.Orien);
        minI = min(output.Infor.Orien);
        
        xgrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
        h{i_off}=line(xgrid, fitS{i_mice,1}{i_off,1}.modelFun(fitS{i_mice,1}{i_off,1}.coefEsts, xgrid), 'Color',colorchoice{i_off});
        hold on;
        plot(fitS{i_mice,1}{i_off,1}.intensityX,fitS{i_mice,1}{i_off,1}.fractCorrY, 'o','Color',colorchoice{i_off});
        %thresh = coefEsts(1)*[1 1];
        plot(fitS{i_mice,1}{i_off,1}.thresh*[1 1], [0 fitS{i_mice,1}{i_off,1}.threshY], '--','Color',colorchoice{i_off});
        plot(bootStats{i_mice,1}{i_off,1}.ci95, fitS{i_mice,1}{i_off,1}.threshY*[1 1], 'Color',colorchoice{i_off});
        
        % set limits correctly
        xLim = [min(xgrid) max(xgrid)].* [0.75 1.25];
        xLim = 10.^ceil(log10(xLim) - [1 0]);
        
        
    end
    legend([h{1} h{2} h{3}], '250ms','500ms','750ms','Location','best')
    legend('boxoff')
    if i_mice==1
        a=5;
    else
        a=10;
    end
    axis([a 100 0 1])
    
    set(gca,'xscale','log','XTick',[10:10:100])
    set(gca,'XTickLabel',{'10','','','', '50','','','','','100'})
    set(gca,'TickDir','out')
    xlabel('Orientation change degree')
    ylabel('Hit Rate')
    title([Allmice{i_mice} 'Hit Rate'])
    
end
%% hit rate seperate by pre offs
figure
h={};
colorchoice = {[0 0 0] [0.5 0.5 0.5] [1 0.8 0.4] [1 0 0]};
for  i_mice = 1:size(Allmice,2)
    subplot(2,2,i_mice)
    for i_off =1:3 %[1,3]
        
        maxI = max(output.Infor.Orien);
        minI = min(output.Infor.Orien);
        
        xgrid = logspace(log10(minI*0.1),log10(maxI*1.5),100);
        h{i_off}=line(xgrid, fitPre{i_mice,1}{i_off,1}.modelFun(fitPre{i_mice,1}{i_off,1}.coefEsts, xgrid), 'Color',colorchoice{i_off});
        
        hold on;
        plot(fitPre{i_mice,1}{i_off,1}.intensityX,fitPre{i_mice,1}{i_off,1}.fractCorrY, 'o','Color',colorchoice{i_off});
        %thresh = coefEsts(1)*[1 1];
        plot(fitPre{i_mice,1}{i_off,1}.thresh*[1 1], [0 fitPre{i_mice,1}{i_off,1}.threshY], '--','Color',colorchoice{i_off});
        plot(bootStatspre{i_mice,1}{i_off,1}.ci95, fitPre{i_mice,1}{i_off,1}.threshY*[1 1], 'Color',colorchoice{i_off});
        
        % set limits correctly
        xLim = [min(xgrid) max(xgrid)].* [0.75 1.25];
        xLim = 10.^ceil(log10(xLim) - [1 0]);
        
        
    end
    legend([h{1} h{2} h{3}], '250ms','500ms','750ms','Location','best')
    legend('boxoff')
    if i_mice==1
        a=5;
    else
        a=10;
    end
    axis([a 93 0 1.05])
    
    set(gca,'xscale','log','XTick',[5 10 20 30 50 90])
    xlabel('Orientation change degree')
    ylabel('Hit Rate')
    title([Allmice{i_mice} 'Hit Rate'])
    
end

%% plot the threshod change seperate by n offs, for all n-1 trials
% Thresh(i_mice,i_off)
% Threshpre(i_mice,i_off)
%
% Thresh_sep(i_mice,i_off,ii_off)



T_ratio = Thresh./repmat(Thresh(:,1),1,3);
T_preratio = Threshpre./repmat(Threshpre(:,1),1,3);


figure
subplot(2,2,1)
for i_mice = 1:size(Allmice,2)
    h=scatter(output.Infor.Off,T_ratio(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:),'MarkerFaceColor',Colors(i_mice,:));
    h.SizeData = 40;
    hold on
    %  h= scatter(output.Infor.Off,T_preratio(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:));
    %  h.SizeData = 55;
    text(220,i_mice*0.1,Allmice{i_mice},'color',Colors(i_mice,:),'FontSize',12)
    ylim([0 1.3])
    xlim([200 800])
end
hold on
h=errorbar(output.Infor.Off,mean(T_ratio,1),std(T_ratio,[],1)./sqrt(size(Allmice,2)),'k');
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'o';
hold on
hline(1,'k:')
ylabel('Orien Thresh Ratio')
set(gca,'XTick',output.Infor.Off)
xlabel('ISI(ms)-n offs')
title('Threshhold-all n-1 trials')

for i_off = 1:3
    subplot(2,2,i_off+1)
    T_ratio_sep=[];
    temp = squeeze(Thresh_sep(:,i_off,:));
    T_ratio_sep = temp./repmat(temp(:,1),1,3);
    for i_mice = 1:size(Allmice,2)
        h=scatter(output.Infor.Off,T_ratio_sep(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:),'MarkerFaceColor',Colors(i_mice,:));
        h.SizeData = 40;
        hold on
        %  h= scatter(output.Infor.Off,T_preratio(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:));
        %  h.SizeData = 55;
        %text(220,i_mice*0.1,Allmice{i_mice},'color',Colors(i_mice,:),'FontSize',12)
        ylim([0 1.3])
        xlim([200 800])
    end
    hold on
    h=errorbar(output.Infor.Off,mean(T_ratio_sep,1),std(T_ratio_sep,[],1)./sqrt(size(Allmice,2)),'k');
    h.LineStyle = 'none';
    h.LineWidth = 1.3;
    h.Marker = 'o';
    hold on
    hline(1,'k:')
    ylabel('Orien Thresh Ratio')
    set(gca,'XTick',output.Infor.Off)
    xlabel('ISI(ms)-n offs')
    title(['Threshhold-' num2str(output.Infor.Off(i_off)) 'n-1 trials'])
    
end
%% plot the threshod change seperate by n-1 offs, for all n trials
% Thresh(i_mice,i_off)
% Threshpre(i_mice,i_off)
%
% Thresh_sep(i_mice,i_off,ii_off)



T_ratio = Thresh./repmat(Thresh(:,1),1,3);
T_preratio = Threshpre./repmat(Threshpre(:,1),1,3);


figure
subplot(2,2,1)
for i_mice = 1:size(Allmice,2)
    h=scatter(output.Infor.Off,T_preratio(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:),'MarkerFaceColor',Colors(i_mice,:));
    h.SizeData = 40;
    hold on
    %  h= scatter(output.Infor.Off,T_preratio(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:));
    %  h.SizeData = 55;
    text(220,i_mice*0.1,Allmice{i_mice},'color',Colors(i_mice,:),'FontSize',12)
    ylim([0 1.3])
    xlim([200 800])
end
hold on
h=errorbar(output.Infor.Off,mean(T_preratio,1),std(T_preratio,[],1)./sqrt(size(Allmice,2)),'k');
h.LineStyle = 'none';
h.LineWidth = 1.3;
h.Marker = 'o';
hold on
hline(1,'k:')
ylabel('Orien Thresh Ratio')
set(gca,'XTick',output.Infor.Off)
xlabel('ISI(ms)-n-1 offs')
title('Threshhold-all n trials')

for i_off = 1:3
    subplot(2,2,i_off+1)
    T_ratio_sep=[];
    temp = squeeze(Thresh_sep(:,:,i_off));
    T_ratio_sep = temp./repmat(temp(:,1),1,3);
    for i_mice = 1:size(Allmice,2)
        h=scatter(output.Infor.Off,T_ratio_sep(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:),'MarkerFaceColor',Colors(i_mice,:));
        h.SizeData = 40;
        hold on
        %  h= scatter(output.Infor.Off,T_preratio(i_mice,:),'MarkerEdgeColor',Colors(i_mice,:));
        %  h.SizeData = 55;
        %text(220,i_mice*0.1,Allmice{i_mice},'color',Colors(i_mice,:),'FontSize',12)
        ylim([0 1.3])
        xlim([200 800])
    end
    hold on
    h=errorbar(output.Infor.Off,mean(T_ratio_sep,1),std(T_ratio_sep,[],1)./sqrt(size(Allmice,2)),'k');
    h.LineStyle = 'none';
    h.LineWidth = 1.3;
    h.Marker = 'o';
    hold on
    hline(1,'k:')
    ylabel('Orien Thresh Ratio')
    set(gca,'XTick',output.Infor.Off)
    xlabel('ISI(ms)-n-1 offs')
    title(['Threshhold-' num2str(output.Infor.Off(i_off)) 'n trials'])
    
end

%% plot normalized thred when fixed N:750ms 
i_off = 3;
T_ratio_sep=[];
    temp = squeeze(Thresh_sep(:,:,i_off));
    T_ratio_sep = temp./repmat(temp(:,1),1,3);
    for i_mice = 1:size(Allmice,2)
        scatter(output.Infor.Off,T_ratio_sep(i_mice,:),'MarkerEdgeColor',[0.5 0.5 0.5],'SizeData',20);
       hold on
        ylim([0 1.4])
        xlim([200 800])
    end
    hold on
    h=errorbar(output.Infor.Off,mean(T_ratio_sep,1),std(T_ratio_sep,[],1)./sqrt(size(Allmice,2)),'k');
    h.LineStyle = 'none';
    h.LineWidth = 1.3;
    h.Marker = 'none';
    hold on
    scatter(output.Infor.Off,mean(T_ratio_sep,1),'MarkerEdgeColor',[0 0 0],'SizeData',60)
    hline(1,'k:')
    ylabel('Orien Thresh Ratio')
    set(gca,'XTick',output.Infor.Off)
    set(gca,'TickDir','out')
    xlabel('ISI(ms)-n-1 offs')
    title(['Threshhold-' num2str(output.Infor.Off(i_off)) 'n trials'])
    %% plot FA rate and threshold on sequential N-1 offs, when collapsed all Ns
   
    figure
    subplot(1,2,1)
    T_preratio = Thresh_prebin./repmat(Thresh_prebin(:,1),1,3);
    for i_mice = 1:size(Allmice,2)
        scatter(output.Infor.Off,T_preratio(i_mice,:),'MarkerEdgeColor',[0.5 0.5 0.5],'SizeData',20);
        hold on
        ylim([0 1.2])
        xlim([200 800])
    end
    hold on
    h=errorbar(output.Infor.Off,mean(T_preratio,1),std(T_preratio,[],1)./sqrt(size(Allmice,2)),'k');
    h.LineStyle = 'none';
    h.LineWidth = 1.3;
    h.Marker = 'none';
    hold on
    scatter(output.Infor.Off,mean(T_preratio,1),'MarkerEdgeColor',[0 0 0],'SizeData',60)
    hline(1,'k:')
    ylabel('Orien Thresh Ratio')
    set(gca,'XTick',output.Infor.Off)
    set(gca,'YTick',[0:0.4:1.2])
    set(gca,'TickDir','out')
    xlabel('ISI(ms)-n-1 offs')
    title(['Threshhold-all N-sequential N-1'])
    
 %  [p1,tbl,stats] = anova1(FA_bin_pre2);
%  b= multcompare(stats);
    
    subplot(1,2,2)
    for i_mice = 1:size(Allmice,2)
        scatter(output.Infor.Off,FA_bin_pre2(i_mice,:),'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Sizedata',20)
        hold on
        
        ylim([0 0.12])
        xlim([200 800])
    end
    hold on
    h=errorbar(output.Infor.Off,mean(FA_bin_pre2,1),std(FA_bin_pre2,[],1)./sqrt(size(Allmice,2)),'k');
    h.LineStyle = 'none';
    h.LineWidth = 1.3;
    h.Marker = 'none';
    
    scatter(output.Infor.Off,mean(FA_bin_pre2,1),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'Sizedata',60)
    
    ylabel('FA rate ')
    set(gca,'XTick',output.Infor.Off)
    set(gca,'YTick',[0:0.03:0.12])
    set(gca,'TickDir','out')
    xlabel('ISI(ms)- n-1 offs')
    title(['FA rate-all N-sequential N-1'])
    
    
    
    
    %% plot FA rate and threshold on sequential N-1 offs, when  fixed N: 750ms 
figure
i=1;
for i_off = 3%1:length(output.Infor.Off)
subplot(1,2,i)
for i_mice = 1:size(Allmice,2)
        scatter(output.Infor.Off,squeeze(Avg_pre_FA(:,i_off,i_mice)),'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'Sizedata',20)
        hold on
       
        ylim([0 0.12])
        xlim([200 800])
    end
    hold on
    h=errorbar(output.Infor.Off,mean(squeeze(Avg_pre_FA(:,i_off,:)),2),std(squeeze(Avg_pre_FA(:,i_off,:)),[],2)./sqrt(size(Allmice,2)),'k');
    h.LineStyle = 'none';
    h.LineWidth = 1.3;
    h.Marker = 'none';
    
    scatter(output.Infor.Off,mean(squeeze(Avg_pre_FA(:,i_off,:)),2),'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1],'Sizedata',60)
    
    ylabel('FA rate ')
    set(gca,'XTick',output.Infor.Off)
    set(gca,'YTick',[0:0.03:0.12])
    set(gca,'TickDir','out')
    xlabel('ISI(ms)- n-1 offs')
    title(['FA rate-' num2str(output.Infor.Off(i_off))  'N' 'Sequential N-1'])
    i=i+1;
   subplot(1,2,i) 
   
     T_ratio_sep=[];
    temp = squeeze(Avg_pre_thresh(:,i_off,:))';
    T_ratio_sep = temp./repmat(temp(:,1),1,3);
    for i_mice = 1:size(Allmice,2)
        scatter(output.Infor.Off,T_ratio_sep(i_mice,:),'MarkerEdgeColor',[0.5 0.5 0.5],'SizeData',20);
       hold on
         ylim([0 max(T_ratio_sep(:))+0.1])
        xlim([200 800])
    end
    hold on
    h=errorbar(output.Infor.Off,mean(T_ratio_sep,1),std(T_ratio_sep,[],1)./sqrt(size(Allmice,2)),'k');
    h.LineStyle = 'none';
    h.LineWidth = 1.3;
    h.Marker = 'none';
    hold on
    scatter(output.Infor.Off,mean(T_ratio_sep,1),'MarkerEdgeColor',[0 0 0],'SizeData',60)
    hline(1,'k:')
    ylabel('Orien Thresh Ratio')
    set(gca,'XTick',output.Infor.Off)
    set(gca,'TickDir','out')
    xlabel('ISI(ms)-n-1 offs')
    title(['Threshhold-' num2str(output.Infor.Off(i_off)) 'N' 'Sequential N-1'])
    i=i+1;
end
%% plot the all trials and matched trials for FA and threshold in the same sale
figure
colors = [0.8 0.8 0.8; 0.5 0.5 0.5;0.3 0.3 0.3];
FAcolors =brewermap(3,'greens');
xmean = mean(Length_FA,1)./1000;
ymean = mean(FA,1);
xerror = std(Length_FA,[],1)./(1000*sqrt(size(Allmice,2)));
yerror = std(FA,[],1)./sqrt(size(Allmice,2));
yyaxis left
for i_off =1:3
errorbar(xmean(i_off),ymean(i_off),yerror(i_off),yerror(i_off),xerror(i_off),xerror(i_off),'color',FAcolors(i_off,:),'LineStyle','none','Marker','o')
hold on
end 
ylim([0 0.12])
set(gca,'YTick',0:0.03:0.12,'TickDir','out')
ylabel('FA rate')

Thresh_N = [];
Thresh_N = Thresh./repmat(Thresh(:,1),1,3);
xmean = mean(Length_Hit,1)./1000;
ymean = mean(Thresh_N,1);
xerror = std(Length_Hit,[],1)./(1000*sqrt(size(Allmice,2)));
yerror = std(Thresh_N,[],1)./sqrt(size(Allmice,2));
yyaxis right
for i_off =1:3
errorbar(xmean(i_off),ymean(i_off),yerror(i_off),yerror(i_off),xerror(i_off),xerror(i_off),'color',colors(i_off,:),'LineStyle','none','Marker','o')
hold on
end 
ylim([0 1])
set(gca,'YTick',0:0.1:1,'TickDir','out')
xlim([1.9 3.6])
set(gca,'XTick',2:0.5:3.5)
ylabel('Threshold')
xlabel('Trial length(s)')

%% get the matched
figure
colors = [0.8 0.8 0.8; 0.5 0.5 0.5;0.3 0.3 0.3];
FAcolors =brewermap(3,'greens');
xmean = mean(MLength_FA,1)./1000;
ymean = mean(MFA,1);
xerror = std(MLength_FA,[],1)./(1000*sqrt(size(Allmice,2)));
yerror = std(MFA,[],1)./sqrt(size(Allmice,2));
yyaxis left
for i_off =1:3
errorbar(xmean(i_off),ymean(i_off),yerror(i_off),yerror(i_off),xerror(i_off),xerror(i_off),'color',FAcolors(i_off,:),'LineStyle','none','Marker','o')
hold on
end 
ylim([0 0.12])
set(gca,'YTick',0:0.03:0.12,'TickDir','out')
ylabel('FA rate')

Thresh_N = [];
Thresh_N = MThresh./repmat(MThresh(:,1),1,3);
xmean = mean(MLength_Hit,1)./1000;
ymean = mean(Thresh_N,1);
xerror = std(MLength_Hit,[],1)./(1000*sqrt(size(Allmice,2)));
yerror = std(Thresh_N,[],1)./sqrt(size(Allmice,2));
yyaxis right
for i_off =1:3
errorbar(xmean(i_off),ymean(i_off),yerror(i_off),yerror(i_off),xerror(i_off),xerror(i_off),'color',colors(i_off,:),'LineStyle','none','Marker','o')
hold on
end 
ylim([0 1])
set(gca,'YTick',0:0.1:1,'TickDir','out')
xlim([1.9 3.6])
set(gca,'XTick',2:0.5:3.5)
ylabel('Threshold')
xlabel('Trial length(s)')

%% plot the all offs and matched offs for normalized threshold
 p = anova1(Length_Hit)
 
figure
subplot(2,2,1)

Thresh_N = [];
Thresh_N = Thresh./repmat(Thresh(:,1),1,3);
xmean = mean(Length_Hit,1)./1000;
ymean = mean(Thresh_N,1);
xerror = std(Length_Hit,[],1)./(1000*sqrt(size(Allmice,2)));
yerror = std(Thresh_N,[],1)./sqrt(size(Allmice,2));
for i_off =1:3
errorbar(xmean(i_off),ymean(i_off),yerror(i_off),yerror(i_off),xerror(i_off),xerror(i_off),'color',colors(i_off,:),'LineStyle','none','Marker','o')
hold on
end 
ylim([0 1.1])
xlim([2.6 3.6])
ylabel('Thresh Norm to 250ms')
set(gca,'XTick',2.6:0.2:3.6)
set(gca,'TickDir','out')
xlabel('Trial length (s)')
title('all trials')

subplot(2,2,2)

Thresh_N = [];
Thresh_N = MThresh./repmat(MThresh(:,1),1,3);
xmean = mean(MLength_Hit,1)./1000;
ymean = mean(Thresh_N,1);
xerror = std(MLength_Hit,[],1)./(1000*sqrt(size(Allmice,2)));
yerror = std(Thresh_N,[],1)./sqrt(size(Allmice,2));
for i_off =1:3
errorbar(xmean(i_off),ymean(i_off),yerror(i_off),yerror(i_off),xerror(i_off),xerror(i_off),'color',colors(i_off,:),'LineStyle','none','Marker','o')
hold on
end
ylim([0 1.1])
xlim([2.6 3.6])
ylabel('Thresh Norm to 250ms')
set(gca,'XTick',2.6:0.2:3.6)
set(gca,'TickDir','out')
xlabel('Trial length (s)')
title('matched trials')
% plot FA rate
subplot(2,2,3)

xmean = mean(Length_FA,1)./1000;
ymean = mean(FA,1);
xerror = std(Length_FA,[],1)./(1000*sqrt(size(Allmice,2)));
yerror = std(FA,[],1)./sqrt(size(Allmice,2));
for i_off =1:3
errorbar(xmean(i_off),ymean(i_off),yerror(i_off),yerror(i_off),xerror(i_off),xerror(i_off),'color',colors(i_off,:),'LineStyle','none','Marker','o')
hold on
end 
ylim([0 0.12])
xlim([1.8 2.8])
ylabel('FA rate')
set(gca,'XTick',1.8:0.2:2.8)
set(gca,'TickDir','out')
xlabel('Trial length (s)')
title('all trials')

% plot FA rate
subplot(2,2,4)

xmean = mean(MLength_FA,1)./1000;
ymean = mean(MFA,1);
xerror = std(MLength_FA,[],1)./(1000*sqrt(size(Allmice,2)));
yerror = std(MFA,[],1)./sqrt(size(Allmice,2));
for i_off =1:3
errorbar(xmean(i_off),ymean(i_off),yerror(i_off),yerror(i_off),xerror(i_off),xerror(i_off),'color',colors(i_off,:),'LineStyle','none','Marker','o')
hold on
end 
ylim([0 0.12])
xlim([1.8 2.8])
ylabel('FA rate')
set(gca,'XTick',1.8:0.2:2.8)
set(gca,'TickDir','out')
xlabel('Trial length (s)')
title('matched trials')
%% 

figure
subplot(2,1,2)
shadedErrorBar_ch(cbin_FA./1000, mean(FA_Trl,1),std(FA_Trl,[],1)./sqrt(size(Allmice,2)))

ylim([0 0.12])
xlim([0 6])
ylabel('FA rate')
set(gca,'XTick',0:2:6)
set(gca,'YTick',0:0.03:0.12,'TickDir','out')
xlabel('Trial length (s)')
subplot(2,1,1)
Thresh_Trl_N =[];
Thresh_Trl_N = Thresh_Trl./repmat(Thresh_Trl(:,1),1,length(cbin_Hit));
Thresh_N = Thresh./repmat(Thresh_Trl(:,1),1,3);


shadedErrorBar_ch(cbin_Hit./1000, mean(Thresh_Trl_N ,1),std(Thresh_Trl_N ,[],1)./sqrt(size(Allmice,2)))

ylim([0 1])
xlim([0 7])
ylabel('Thresh Norm to shortest length')
set(gca,'XTick',0:2:6)
set(gca,'TickDir','out')
xlabel('Trial length (s)')

%% plot the trial length dependence of FA rate and Thresh
figure
subplot(2,1,1)
shadedErrorBar_ch(cbin_FA./1000, mean(FA_Trl,1),std(FA_Trl,[],1)./sqrt(size(Allmice,2)))
hold on
scatter(mean(Length_FA,1)./1000,mean(FA,1),'k')
errorbar(mean(Length_FA,1)./1000,mean(FA,1), std(FA,[],1)./sqrt(size(Allmice,2)),'k')
ylim([0 0.12])
ylabel('FA rate')
set(gca,'XTick',0:0.5:6)
set(gca,'TickDir','out')
xlabel('Trial length (s)')
subplot(2,1,2)
Thresh_Trl_N =[];
Thresh_Trl_N = Thresh_Trl./repmat(Thresh_Trl(:,1),1,length(cbin_Hit));
Thresh_N = Thresh./repmat(Thresh_Trl(:,1),1,3);


shadedErrorBar_ch(cbin_Hit./1000, mean(Thresh_Trl_N ,1),std(Thresh_Trl_N ,[],1)./sqrt(size(Allmice,2)))
hold on
scatter(mean(Length_Hit,1)./1000,mean(Thresh_N,1),'k')
errorbar(mean(Length_Hit,1)./1000,mean(Thresh_N,1), std(Thresh_N,[],1)./sqrt(size(Allmice,2)),'k')
ylim([0 1])
ylabel('Thresh Norm to shortest length')
set(gca,'XTick',0:0.5:7)
set(gca,'TickDir','out')
xlabel('Trial length (s)')



p = anova1(Thresh_Trl_N)
%% plot c, dprime, hit 

colors = {[0.2 0.6 1], [0 0 0],[1 0 0]};
figure
subplot(2,2,1)
temp = [];
temp = cat(3,FA-repmat(FA(:,2),1,3),hit_22Deg -repmat(hit_22Deg(:,2),1,3), hit_90Deg -repmat(hit_90Deg(:,2),1,3));

for i_off = 1:3
    errorbar([0 22.5 90],mean(squeeze(temp(:,i_off,:)),1), std(squeeze(temp(:,i_off,:)),[],1)./sqrt(size(Allmice,2)),'Marker','o','Color',colors{i_off});
    hold on
end 

xlim([-5 95])
ylim([-0.15 0.15])
hline(0,'k:')
set(gca,'XTick',0:45:90,'YTick',-0.15:0.05:0.15,'TickDir','out')
xlabel('Orientation (deg)')
ylabel('\Delta Hit')
axis square

subplot(2,2,2)

temp = [];
temp = cat(3,c_22Deg -repmat(c_22Deg(:,2),1,3), c_90Deg -repmat(c_90Deg(:,2),1,3));

for i_off = 1:3
    errorbar([22.5 90],mean(squeeze(temp(:,i_off,:)),1), std(squeeze(temp(:,i_off,:)),[],1)./sqrt(size(Allmice,2)),'Marker','o','Color',colors{i_off});
    hold on
end 

xlim([-5 95])
ylim([-0.4 0.4])
hline(0,'k:')
set(gca,'XTick',0:45:90,'YTick',-0.4:0.2:0.4,'TickDir','out')
xlabel('Orientation (deg)')
ylabel('\Delta c')
axis square

subplot(2,2,3)

temp = [];
temp = cat(3,d_22Deg -repmat(d_22Deg(:,2),1,3), d_90Deg -repmat(d_90Deg(:,2),1,3));

for i_off = 1:3
    errorbar([22.5 90],mean(squeeze(temp(:,i_off,:)),1), std(squeeze(temp(:,i_off,:)),[],1)./sqrt(size(Allmice,2)),'Marker','o','Color',colors{i_off});
    hold on
end 

xlim([-5 95])
ylim([-0.4 0.4])
hline(0,'k:')
set(gca,'XTick',0:45:90,'YTick',-0.4:0.2:0.4,'TickDir','out')
xlabel('Orientation (deg)')
ylabel('\Delta d prime')
axis square
%% plot only c changes 
temp = [];
temp = cat(3,c_22Deg -repmat(c_22Deg(:,2),1,3), c_90Deg -repmat(c_90Deg(:,2),1,3));
d_c = squeeze(temp(:,:,1));
Num = size(temp,1);
figure
scatter(repmat(0,1,Num),d_c(:,1),'SizeData',10, 'MarkerEdgeColor',[0.5 0.5 0.5])
hold on
errorbar(0,mean(d_c(:,1)),std(d_c(:,1))./sqrt(Num),'Color',[0.2 0.6 1],'Marker','o')

scatter(repmat(0.5,1,Num),d_c(:,3),'SizeData',10, 'MarkerEdgeColor',[0.5 0.5 0.5])

errorbar(0.5,mean(d_c(:,3)),std(d_c(:,3))./sqrt(Num),'Color',[1 0 0],'Marker','o')

xlim([-0.3 0.8])
ylim([-0.5 0.5])

hline(0,'k:')
set(gca,'YTick',-0.5:0.25:0.5,'XTick',[0 0.5],'TickDir','out')
ylabel('\Delta criterion')
axis square

