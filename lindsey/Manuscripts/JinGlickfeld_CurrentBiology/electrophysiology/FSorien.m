%% new pair pulse analysis
function  [FS]=FSorien(input,Trigger,spikes)
%%
global Fs stimOnTimeS Offs
Fs=double(Trigger.SamplingFreq);

[Time_start_new,Time_stop_new,photo_duration_MS,photo_inter_MS,~]=PhotoISI(Trigger,input);
%[Time_start_new,Time_stop_new,photo_duration_MS,photo_inter_MS,~]=PhotoISI_nolabjack(Trigger,input);
% use the LFP identify the layers of each template
baseline=0.2;
binwidth = 0.020;

[Data,psth,~,~,psth_full,spon]=align_isi(spikes,baseline,binwidth,Time_start_new); % add peak+-50ms area for stat test 



%% now separate by different off conditions of the psth_full
StimOffs = double(cell2mat(input.tStimOffTimeMs));
Offs = unique(StimOffs);
edges = -baseline:binwidth:0.35;
base = find(edges==0); % use this to subtract the peak response
basebins = round(baseline./binwidth);
Ix_window = edges>=0 & edges <= 0.13; % for the timewindow 120 ms during the visual presentation
on_edges =  edges(Ix_window);
cycle =  unique(cell2mat(input.tCyclesOn))+1;
%trialnum = length(psth.First{1,1});
B_Orien =double( cell2mat(input.tBaseGratingDirectionDeg));
%B_Orien = B_Orien(1:trialnum);
BOrien = unique(B_Orien);

T_Orien = cell2mat(input.tGratingDirectionDeg);
T_Orien = round(mod(T_Orien+180,180).*10,0)./10;
%T_Orien = T_Orien(1:trialnum);
TOrien = unique(T_Orien);
FS.BOrien = BOrien;
FS.theta = deg2rad(BOrien);
FS.fitorien = 0:1:179;
FS.fittheta = deg2rad(FS.fitorien);
% concate of the two  offs
for i_template = 1:length(psth.area)
    
    % find the peak of the First,target and adpt all apply it to all the pulses
    
    
    First = [];
    First = cell2mat(cellfun(@(x) x(1,:),psth.raw{i_template,1},'UniformOutput',false));
    Target = [];
    Target = cell2mat(cellfun(@(x) x(6,:),psth.raw{i_template,1},'UniformOutput',false));
    Adpt_4 = [];
    Adpt_4 = cell2mat(cellfun(@(x) x(4,:),psth.raw{i_template,1},'UniformOutput',false));
    Adpt_5 = [];
    Adpt_5 = cell2mat(cellfun(@(x) x(5,:),psth.raw{i_template,1},'UniformOutput',false));
    
    [~, ~,FS.responseIX(i_template,1)]=findpeakmiao((mean(First,1)+mean(Target,1)+mean([Adpt_4;Adpt_5],1))./3,Ix_window,on_edges,basebins);
    
    for i_orien = 1:length(BOrien)
        % if it is significant driven by any of the orientations
        [FS.stats_ori(i_template,i_orien),FS.stats_oriP(i_template,i_orien)]=ttest2(psth.First{i_template}(B_Orien==BOrien(i_orien)),psth.baseline{i_template}(B_Orien==BOrien(i_orien)),'Tail','right'); %fixed window
        [FS.stats_oriN(i_template,i_orien),FS.stats_oriPN(i_template,i_orien)]=ttest2(psth.areaN{i_template,1}(B_Orien==BOrien(i_orien),1),psth.baseline{i_template}(B_Orien==BOrien(i_orien)),'Tail','right'); %adjusted window
        FS.Traces.first.raw_mean(i_template,i_orien,:) = mean(First(B_Orien==BOrien(i_orien),:),1); % first response traces
        FS.Traces.first.raw_sem(i_template,:) = std(First(B_Orien==BOrien(i_orien),:),[],1)./sqrt(sum(B_Orien==BOrien(i_orien)));
        
        % organize in the same orientation order
        
        FS.Traces.target.raw_mean(i_template,i_orien,:) = mean(Target(T_Orien==BOrien(i_orien),:),1); % first response traces
        FS.Traces.target.raw_sem(i_template,i_orien,:) = std(Target(T_Orien==BOrien(i_orien),:),[],1)./sqrt(sum(T_Orien==BOrien(i_orien)));
        
        
        FS.Traces.adpt.raw_mean(i_template,i_orien,:) = mean([Adpt_4(B_Orien==BOrien(i_orien),:); Adpt_5(B_Orien==BOrien(i_orien),:) ],1); % first response traces
        FS.Traces.adpt.raw_sem(i_template,i_orien,:) = std([Adpt_4(B_Orien==BOrien(i_orien),:); Adpt_5(B_Orien==BOrien(i_orien),:) ],[],1)./sqrt(sum(B_Orien==BOrien(i_orien))*2);
        
        % get the response raw that subtracted by bin 0
        FS.response.first.raw{i_template,i_orien} = First(B_Orien==BOrien(i_orien),FS.responseIX(i_template,1))-mean(First(B_Orien==BOrien(i_orien),base)); % already subtracted by bin 0
        FS.response.first.mean(i_template,i_orien) = mean(First(B_Orien==BOrien(i_orien),FS.responseIX(i_template,1))-mean(First(B_Orien==BOrien(i_orien),base)));
        FS.response.first.sem(i_template,i_orien) = std(First(B_Orien==BOrien(i_orien),FS.responseIX(i_template,1))-mean(First(B_Orien==BOrien(i_orien),base)))./sqrt(sum(B_Orien==BOrien(i_orien)));
        % get the area response for both adapt and first
        
        FS.response.first.area_raw{i_template,i_orien} = psth.areaN{i_template,1}(B_Orien==BOrien(i_orien),1)-mean(First(B_Orien==BOrien(i_orien),base));
        FS.response.first.area_mean(i_template,i_orien) = mean(FS.response.first.area_raw{i_template,i_orien});
        FS.response.first.area_sem(i_template,i_orien) = std(FS.response.first.area_raw{i_template,i_orien})./sqrt(sum(B_Orien==BOrien(i_orien)));
        
        
        FS.response.adpt.area_raw{i_template,i_orien} = psth.areaN{i_template,1}(B_Orien==BOrien(i_orien),2)-mean(First(B_Orien==BOrien(i_orien),base));
        FS.response.adpt.area_mean(i_template,i_orien) = mean(FS.response.adpt.area_raw{i_template,i_orien});
        FS.response.adpt.area_sem(i_template,i_orien) = std(FS.response.adpt.area_raw{i_template,i_orien})./sqrt(sum(B_Orien==BOrien(i_orien)));
        
        
        % get the adapt response
        temp=[];
        idx_temp = [];
        temp = [Adpt_4; Adpt_5];
        idx_temp = [B_Orien==BOrien(i_orien) B_Orien==BOrien(i_orien)];
        FS.response.adpt.raw{i_template,i_orien} =  temp(idx_temp,FS.responseIX(i_template,1)) - mean(temp(idx_temp,base));
        Temp=[];
        Temp =  temp(idx_temp,FS.responseIX(i_template,1)) - mean(temp(idx_temp,base));
        FS.response.adpt.meanraw{i_template,i_orien} = (Temp(1:(length(Temp)./2))+Temp(((length(Temp)./2)+1):end))./2; 
        FS.response.adpt.mean(i_template,i_orien) = mean(temp(idx_temp,FS.responseIX(i_template,1)) - mean(temp(idx_temp,base)));
        FS.response.adpt.sem(i_template,i_orien) = std(temp(idx_temp,FS.responseIX(i_template,1)) - mean(temp(idx_temp,base)))./sqrt(sum(idx_temp));
        
        
        % get the target response
        FS.response.target.raw{i_template,i_orien} = Target(T_Orien==BOrien(i_orien),FS.responseIX(i_template,1))-mean(Target(T_Orien==BOrien(i_orien),base)); % already subtracted by bin 0
        FS.response.target.mean(i_template,i_orien) = mean(Target(T_Orien==BOrien(i_orien),FS.responseIX(i_template,1))-mean(Target(T_Orien==BOrien(i_orien),base)));
        FS.response.target.sem(i_template,i_orien) = std(Target(T_Orien==BOrien(i_orien),FS.responseIX(i_template,1))-mean(Target(T_Orien==BOrien(i_orien),base)))./sqrt(sum(T_Orien==BOrien(i_orien)));
        
        % get the adaptation index for each orientation
        FS.adpt_idx(i_template,i_orien) = (FS.response.adpt.mean(i_template,i_orien)./abs(FS.response.first.mean(i_template,i_orien)))-1; 
        FS.target_adpt(i_template,i_orien) =  (FS.response.target.mean(i_template,i_orien)-FS.response.adpt.mean(i_template,i_orien))./abs(FS.response.adpt.mean(i_template,i_orien)); 
        FS.target_first(i_template,i_orien) = (FS.response.target.mean(i_template,i_orien)-FS.response.first.mean(i_template,i_orien))./abs(FS.response.first.mean(i_template,i_orien)); 
        FS.area_adpt_idx(i_template,i_orien) = (FS.response.adpt.area_mean(i_template,i_orien)./abs(FS.response.first.area_mean(i_template,i_orien)))-1; 
        
    end
    
    
    
    for i = 1:1000 % random sample with replacement
        for i_orien = 1:length(BOrien)
            temp = [];
            temp = FS.response.first.raw{i_template,i_orien};
            FS.bootstrap.first(i_template,i,i_orien) =  mean(randsample(temp,length(temp),1));
            temp = [];
            temp = FS.response.adpt.raw{i_template,i_orien};
            FS.bootstrap.adpt(i_template,i,i_orien) =  mean(randsample(temp,length(temp),1));
            temp = [];
            temp = FS.response.target.raw{i_template,i_orien};
            FS.bootstrap.target(i_template,i,i_orien) =  mean(randsample(temp,length(temp),1));
            
            
            
        end
        
        % get the fit for 1000 times
        
    end
    
    % get the data organized as cellsxtrials for glm purpuse,needs to be
    % zscored
    Temp = [];
    Temp = [FS.response.adpt.meanraw{i_template,1}; FS.response.target.raw{i_template,4}];
    Z_temp = [];
    Z_temp = zscore(Temp); 
    FS.glm.resp(:,i_template) = Z_temp;
    
end 
FS.glm.TrlOut = [zeros(1,length(Z_temp)./2) ones(1,length(Z_temp)./2)]; 

for i_template = 1:length(psth.area)
    % determine cells that are significantly driven by 0 deg
   
    [FS.stats_F(i_template),FS.stats_FP(i_template)]=ttest2(psth.First{i_template},psth.baseline{i_template},'Tail','right'); % the mean of X is greater than the mean of Y;
    [FS.stats_areaF(i_template),FS.stats_areaFP(i_template)]=ttest2(psth.areaN{i_template}(:,1),psth.baseline{i_template},'Tail','right'); % the mean of X is greater than the mean of Y;

    [FS.stats_T(i_template),~]=ttest2(psth.target{i_template},psth.target_baseline{i_template},'Tail','right');
    [FS.stats_A(i_template),~]=ttest2([psth.secondbase{i_template,4}; psth.secondbase{i_template,5}],[psth.second{i_template,4}; psth.second{i_template,5}],'Tail','right');
    %
    if FS.stats_F(i_template)==1
        [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(FS.theta,FS.response.first.mean(i_template,:));
        FS.fit.first.R_square(i_template) = R_square;
        if rad2deg(u1_hat)<0
            FS.fit.first.prefer(i_template) = 180 + rad2deg(u1_hat);
        else
            FS.fit.first.prefer(i_template) = rad2deg(u1_hat);
        end
        FS.fit.first.curve(i_template,1:180) = b_hat + R1_hat.*exp(k1_hat.*(cos(2.*(FS.fittheta-u1_hat))-1));
        % calculate OSI, prefer-nonprefer/(prefer+nonprefer)
        p_temp=[];
        o_temp=[];
        p_temp =  b_hat + R1_hat.*exp(k1_hat.*(cos(2.*(u1_hat-u1_hat))-1));
        o_temp =  b_hat + R1_hat.*exp(k1_hat.*(cos(2.*(u1_hat+deg2rad(90)-u1_hat))-1));
        FS.fit.first.OSI(i_template) = (p_temp-o_temp)./(p_temp+o_temp); 
        % integral under the curve
        fun = @(x)  b_hat + R1_hat.*exp(k1_hat.*(cos(2.*(x-u1_hat))-1)); 
        FS.fit.first.integral(i_template) = integral(fun,FS.fittheta(1),FS.fittheta(end));  
        
        
        [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(FS.theta,FS.response.adpt.mean(i_template,:));
        FS.fit.adpt.R_square(i_template) = R_square;
        if rad2deg(u1_hat)<0
            FS.fit.adpt.prefer(i_template) = 180 + rad2deg(u1_hat);
        else
            FS.fit.adpt.prefer(i_template) = rad2deg(u1_hat);
        end
        FS.fit.adpt.curve(i_template,1:180) = b_hat + R1_hat.*exp(k1_hat.*(cos(2.*(FS.fittheta-u1_hat))-1));
        % calculate OSI, prefer-nonprefer/(prefer+nonprefer)
        p_temp=[];
        o_temp=[];
        p_temp =  b_hat + R1_hat.*exp(k1_hat.*(cos(2.*(u1_hat-u1_hat))-1));
        o_temp =  b_hat + R1_hat.*exp(k1_hat.*(cos(2.*(u1_hat+deg2rad(90)-u1_hat))-1));
        FS.fit.adpt.OSI(i_template) = (p_temp-o_temp)./(p_temp+o_temp); 
        % integral under the curve
        fun = @(x)  b_hat + R1_hat.*exp(k1_hat.*(cos(2.*(x-u1_hat))-1)); 
        FS.fit.adpt.integral(i_template) = integral(fun,FS.fittheta(1),FS.fittheta(end));  
        
        
        [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(FS.theta,FS.response.target.mean(i_template,:));
        FS.fit.target.R_square(i_template) = R_square;
        if rad2deg(u1_hat)<0
            FS.fit.target.prefer(i_template) = 180 + rad2deg(u1_hat);
        else
            FS.fit.target.prefer(i_template) = rad2deg(u1_hat);
        end
        FS.fit.target.curve(i_template,1:180) = b_hat + R1_hat.*exp(k1_hat.*(cos(2.*(FS.fittheta-u1_hat))-1));
         % calculate OSI, prefer-nonprefer/(prefer+nonprefer)
        p_temp=[];
        o_temp=[];
        p_temp =  b_hat + R1_hat.*exp(k1_hat.*(cos(2.*(u1_hat-u1_hat))-1));
        o_temp =  b_hat + R1_hat.*exp(k1_hat.*(cos(2.*(u1_hat+deg2rad(90)-u1_hat))-1));
        FS.fit.target.OSI(i_template) = (p_temp-o_temp)./(p_temp+o_temp); 
        % integral under the curve
        fun = @(x)  b_hat + R1_hat.*exp(k1_hat.*(cos(2.*(x-u1_hat))-1)); 
        FS.fit.target.integral(i_template) = integral(fun,FS.fittheta(1),FS.fittheta(end)); 
        
        % determine the adaptation index of the integral
        FS.integral.adpt_idx(i_template) = (FS.fit.adpt.integral(i_template)-FS.fit.first.integral(i_template))./abs(FS.fit.first.integral(i_template));
        FS.integral.target_adpt(i_template) = (FS.fit.target.integral(i_template)-FS.fit.adpt.integral(i_template))./abs(FS.fit.adpt.integral(i_template));
        FS.integral.target_first(i_template) = (FS.fit.target.integral(i_template)-FS.fit.first.integral(i_template))./abs(FS.fit.first.integral(i_template));
        
        % determine if it is reliable fit
        
        for i = 1:1000 
            [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(FS.theta, squeeze(FS.bootstrap.first(i_template,i,:)));
            FS.fit.first.bootstrap.pref(i_template,i) = rad2deg(u1_hat);
            [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(FS.theta, squeeze(FS.bootstrap.adpt(i_template,i,:)));
            FS.fit.adpt.bootstrap.pref(i_template,i) = rad2deg(u1_hat);
            [b_hat, k1_hat, R1_hat,u1_hat,sse,R_square] = miaovonmisesfit_ori(FS.theta, squeeze(FS.bootstrap.target(i_template,i,:)));
            FS.fit.target.bootstrap.pref(i_template,i) = rad2deg(u1_hat);
            
            
        end
        temp = [];
        temp = abs(FS.fit.first.bootstrap.pref(i_template,:)-FS.fit.first.prefer(i_template));
        temp(temp>90) = 180 - temp(temp>90);
        
        FS.fit.first.bootstrap.criterion(i_template) = sum(temp<=30)./1000;
        
        temp = [];
        temp = abs(FS.fit.adpt.bootstrap.pref(i_template,:)-FS.fit.adpt.prefer(i_template));
        temp(temp>90) = 180 - temp(temp>90);
        
        FS.fit.adpt.bootstrap.criterion(i_template) = sum(temp<=30)./1000;
        
        temp = [];
        temp = abs(FS.fit.target.bootstrap.pref(i_template,:)-FS.fit.target.prefer(i_template));
        temp(temp>90) = 180 - temp(temp>90);
        
        FS.fit.target.bootstrap.criterion(i_template) = sum(temp<=30)./1000;
        
        
    else
        FS.fit.first.R_square(i_template)=NaN;
        FS.fit.first.prefer(i_template) = NaN;
        FS.fit.first.curve(i_template,1) = NaN;
        FS.fit.adpt.R_square(i_template)=NaN;
        FS.fit.adpt.prefer(i_template) = NaN;
        FS.fit.adpt.curve(i_template,1) = NaN;
        FS.fit.target.R_square(i_template)=NaN;
        FS.fit.target.prefer(i_template) = NaN;
        FS.fit.target.curve(i_template,1) = NaN;
        
        FS.fit.first.bootstrap.pref(i_template,1) = NaN;
        FS.fit.first.bootstrap.criterion(i_template) = NaN;
        FS.fit.adpt.bootstrap.pref(i_template,1) = NaN;
        FS.fit.adpt.bootstrap.criterion(i_template) = NaN;
        FS.fit.target.bootstrap.pref(i_template,1) = NaN;
        FS.fit.target.bootstrap.criterion(i_template) = NaN;
        
    end
    
    
end
% get the fit for the orientation tuning for the first response and target 
FS.stats_both = FS.stats_F&FS.stats_T; 
FS.stats_Ffit = FS.fit.first.bootstrap.criterion>=0.9; 
FS.stats_Afit = FS.fit.adpt.bootstrap.criterion>=0.9; 
FS.stats_Tfit = FS.fit.target.bootstrap.criterion>=0.9; 

%% linear regreassion fit
for i_template = 1:length(psth.area)
    
    tempX = [];
    tempX = FS.response.first.mean(i_template,:);
    tempY = [];
    tempY = FS.response.adpt.mean(i_template,:);
    tempZ = [];
    tempZ =  FS.response.target.mean(i_template,:);
    
    [R P]= corrcoef(tempX,tempY);   
   
    FS.response.adpt.r(i_template) = R(2,1); % correlation coefficient
    FS.response.adpt.p(i_template) = P(2,1); 
    figure
    subplot(1,2,1)
    scatter(tempX, tempY,'k')
    h=lsline; 
    [p] = polyfit(get(h,'xdata'), get(h,'ydata'),1); 
  
    hold on
    plot([0 1],[0 1],'r')
   
    axis square 
    
    FS.response.adpt.slope(i_template) = p(1);
    FS.response.adpt.intercept(i_template) = p(2); 
    
    % for target response
    [R P]= corrcoef(tempX,tempZ);   
   
    FS.response.target.r(i_template) = R(2,1); % correlation coefficient
    FS.response.target.p(i_template) = P(2,1); 
    subplot(1,2,2)
    scatter(tempX, tempZ,'k')
    h=lsline; 
    [p] = polyfit(get(h,'xdata'), get(h,'ydata'),1); 
  
    hold on
    plot([0 1],[0 1],'r')
   
    axis square 
    
    FS.response.target.slope(i_template) = p(1);
    FS.response.target.intercept(i_template) = p(2); 
    
end 
close all

%% plot to check
idx = find((FS.stats_Ffit& FS.stats_F&FS.stats_Afit&FS.stats_Tfit)==1);
FS.allfit.stats = idx; % reliable fit for first, adpt and target
for i =1: length(idx)
    i_template = idx(i);
        figure
        subplot(2,2,1)
        errorbar(BOrien,FS.response.first.mean(i_template,:)./binwidth,FS.response.first.sem(i_template,:)./binwidth,'k','Marker','o','LineStyle','none')
        hold on
        plot(FS.fitorien,FS.fit.first.curve(i_template,:)./binwidth,'k')
        
        errorbar(BOrien,FS.response.adpt.mean(i_template,:)./binwidth,FS.response.adpt.sem(i_template,:)./binwidth,'Color',[0.2 0.6 0.4],'Marker','o','LineStyle','none')
        plot(FS.fitorien,FS.fit.adpt.curve(i_template,:)./binwidth,'color',[0.2 0.6 0.4])
        
        errorbar(BOrien,FS.response.target.mean(i_template,:)./binwidth,FS.response.target.sem(i_template,:)./binwidth,'Color',[1 0 0],'Marker','o','LineStyle','none')
        plot(FS.fitorien,FS.fit.target.curve(i_template,:)./binwidth,'color',[1 0 0])
        xlim([-10 190])
        set(gca,'XTick',0:30:180,'TickDir','out')
        xlabel('Orien (Deg)')
        ylabel('FR (Hz)')
     % align to the preferred orientation in the center and normalized it
        subplot(2,2,2)
     [R_m,pref]=max(FS.fit.first.curve(i_template,:)); 
     FS.allfit.fitF(i_template,:) = circshift(FS.fit.first.curve(i_template,:),-(pref-90)); % so that preferred orientation is in the center
     FS.allfit.fitA(i_template,:) = circshift(FS.fit.adpt.curve(i_template,:),-(pref-90));
     FS.allfit.fitT(i_template,:) = circshift(FS.fit.target.curve(i_template,:),-(pref-90));
     FS.allfit.BOrien(i_template,:) = BOrien-(pref-90);
     FS.allfit.BOrien(i_template,FS.allfit.BOrien(i_template,:)>180) = FS.allfit.BOrien(i_template,FS.allfit.BOrien(i_template,:)>180)-180;
     temp_idx = [];
     temp_idx = FS.allfit.BOrien(i_template,:)<0;
     
     FS.allfit.BOrien(i_template,temp_idx) = FS.allfit.BOrien(i_template,temp_idx)+180; 
     scatter( FS.allfit.BOrien(i_template,:),FS.response.first.mean(i_template,:)./binwidth,'k')
     hold on
     plot(FS.fitorien,FS.allfit.fitF(i_template,:)./binwidth ,'k')
     scatter( FS.allfit.BOrien(i_template,:),FS.response.adpt.mean(i_template,:)./binwidth,'MarkerEdgeColor',[0.2 0.6 0.4])
     plot(FS.fitorien,FS.allfit.fitA(i_template,:)./binwidth,'color',[0.2 0.6 0.4])
     scatter( FS.allfit.BOrien(i_template,:),FS.response.target.mean(i_template,:)./binwidth,'MarkerEdgeColor',[1 0 0])
     plot(FS.fitorien,FS.allfit.fitT(i_template,:)./binwidth,'color',[1 0 0])
     xlim([-10 190])
     set(gca,'XTick',0:30:180,'TickDir','out')
     xlabel('Orien (Deg)')
     ylabel('FR (Hz)')
     % normalize response
     subplot(2,2,3)
     FS.allfit.norm_fitF(i_template,:) = FS.allfit.fitF(i_template,:)./R_m; 
     FS.allfit.norm_fitA(i_template,:) = FS.allfit.fitA(i_template,:)./R_m; 
     FS.allfit.norm_fitT(i_template,:) = FS.allfit.fitT(i_template,:)./R_m; 
     
     plot(FS.fitorien,FS.allfit.norm_fitF(i_template,:) ,'k')
     hold on
     plot(FS.fitorien,FS.allfit.norm_fitA(i_template,:),'color',[0.2 0.6 0.4])
     plot(FS.fitorien,FS.allfit.norm_fitT(i_template,:),'color',[1 0 0])
     xlim([-10 190])
     set(gca,'XTick',0:30:180,'TickDir','out')
     xlabel('Orien (Deg)')
     ylabel('FR (Hz)')
end
  
    

%% get latency measures
for i_template = 1:length(psth.area)
    
    for i_trial = 1:length(StimOffs)
        for i_cycle=1:6
            if ~isempty ( Data{i_template,1}{i_trial,1}{i_cycle})
                if size(Data{i_template,1}{i_trial,1}{i_cycle},1)==1
                    FS.spikes{i_template,i_cycle}{i_trial,1} = Data{i_template,1}{i_trial,1}{i_cycle};
                else
                    FS.spikes{i_template,i_cycle}{i_trial,1} = Data{i_template,1}{i_trial,1}{i_cycle}';
                end
            else
                FS.spikes{i_template,i_cycle}{i_trial,1} = double.empty(1,0);
            end
        end
    end
end
% getkernel density fit to get the latency measurement
FS.kernel_fix.t = linspace(-0.2,0.35,1000);
for i_template = 1:length(psth.area)
    % stored in FS.pp
    temp = [];
    temp = cat(2,FS.spikes{i_template,1}{:}); % all fisrt response
    % using the fixed kernel bandwith for every unit
    FS.kernel.density(i_template,:) = sskernel(temp,FS.kernel_fix.t,0.01);
    % get the latency using the 0.5 height method and 0.2 height
    Base = []; % the same base apply to the adapt condition
    Max = [];
    Thresh = [];
    Base = mean(FS.kernel.density(i_template,FS.kernel_fix.t>=-0.2 & FS.kernel_fix.t<0)); % the baseline before the visual onset (0.2s)
    if isnan(Base)
        Base = 0;
    end
    % get the rebaselined density
    FS.kernel.shape(i_template,:) = FS.kernel.density(i_template,:)-Base;
    Max = max(FS.kernel.density(i_template,FS.kernel_fix.t>=0 & FS.kernel_fix.t<=0.13)); % peak window as 0-0.13 s
    % when SDF exceeds the level given by base+0.5*(max-base);
    Thresh = Base+0.5*(Max - Base);
    T_chunk = [];
    T_chunk = FS.kernel_fix.t(FS.kernel_fix.t>=0 & FS.kernel_fix.t<=0.13);
    TTemp = [];
    TTemp = FS.kernel.density(i_template,FS.kernel_fix.t>=0 & FS.kernel_fix.t<=0.13)>=Thresh;
    if sum( TTemp)>=1
        Latency = T_chunk(TTemp);
        FS.Latency(i_template) = Latency(1);
    else
        FS.Latency(i_template) = NaN;
    end
    % 20% height
    Thresh = [];
    Thresh = Base+0.2*(Max - Base);
    T_chunk = [];
    T_chunk = FS.kernel_fix.t(FS.kernel_fix.t>=0 & FS.kernel_fix.t<=0.13);
    TTemp = [];
    TTemp = FS.kernel.density(i_template,FS.kernel_fix.t>=0 & FS.kernel_fix.t<=0.13)>=Thresh;
    if sum( TTemp)>=1
        Latency = T_chunk(TTemp);
        FS.Latency_20(i_template) = Latency(1);
    else
        FS.Latency_20(i_template) = NaN;
    end
    
    
end

%% get the layers of template and rpv of each cells, redo layers by mean waveform?
% % for V1 poly
% for i_template = 1:length(psth.area)
%     
%     if  spikes.waveform.pos(i_template,2)<=LFP.layers(1)
%         FS.layers(i_template) = 1;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(1) && spikes.waveform.pos(i_template,2)<=LFP.layers(2)
%         FS.layers(i_template) = 3;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(2) && spikes.waveform.pos(i_template,2)<=LFP.layers(3)
%         FS.layers(i_template) = 4;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(3) && spikes.waveform.pos(i_template,2)<=LFP.layers(4)
%         FS.layers(i_template) = 5;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(4)
%         FS.layers(i_template) = 6;
%     end
%     
%     
%     
% end

% for HVAs, 1: above the initial big sink
% 2: in the initial sink
% 3: in between the superficial and deep sink
% 4: the deep sink or below it
% for SC, 1: above the SGS
% 2:IN SGS / 3:in SO / 4: below SO
% % for LGN, 1: above the sink, 2:in the sink, 3:in the source, 4: below the
% % % source
% for i_template = 1:length(psth.area)
% 
%     if  spikes.waveform.pos(i_template,2)<=LFP.layers(1)
%         FS.layers(i_template) = 1;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(1) && spikes.waveform.pos(i_template,2)<=LFP.layers(2)
%         FS.layers(i_template) = 2;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(2) && spikes.waveform.pos(i_template,2)<=LFP.layers(3)
%         FS.layers(i_template) = 3;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(3)
%         FS.layers(i_template) = 4;
% 
%     end
% 
% 
% 
% end


% for layers accross shanks for HVAs
% for i_template = 1:length(psth.area)
%     i_shank = [];
%     if round(spikes.waveform.pos(i_template,1))<=50 % 0
%         i_shank =1;
%     elseif  round(spikes.waveform.pos(i_template,1)) >=300 && round(spikes.waveform.pos(i_template,1))<=500 % 400
%         i_shank =2;
%     elseif round(spikes.waveform.pos(i_template,1)) >=700 && round(spikes.waveform.pos(i_template,1))<=900 % 800
%         i_shank = 3;
%     elseif round(spikes.waveform.pos(i_template,1)) >=1100 && round(spikes.waveform.pos(i_template,1))<=1300 % 1200
%         i_shank = 4;
%     end
%        
%        %temp(i_template) = i_shank; 
%     if  spikes.waveform.pos(i_template,2)<=LFP.layers(i_shank,1)
%         FS.layers(i_template) = 1;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(i_shank,1) && spikes.waveform.pos(i_template,2)<=LFP.layers(i_shank,2)
%         FS.layers(i_template) = 2;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(i_shank,2) && spikes.waveform.pos(i_template,2)<=LFP.layers(i_shank,3)
%         FS.layers(i_template) = 3;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(i_shank,3)
%         FS.layers(i_template) = 4;
% 
%     end
% 
% 
% 
% end

% for layers accross shanks for V1
% for i_template = 1:length(psth.area)
%     
%     i_shank = [];
%     if round(spikes.waveform.pos(i_template,1))==0
%         i_shank =1;
%     elseif  round(spikes.waveform.pos(i_template,1)) ==400
%         i_shank =2;
%     elseif round(spikes.waveform.pos(i_template,1)) ==800
%         i_shank = 3;
%     elseif round(spikes.waveform.pos(i_template,1)) ==1200
%         i_shank = 4;
%     end
%     
%     if  spikes.waveform.pos(i_template,2)<=LFP.layers(i_shank,1)
%         FS.layers(i_template) = 1;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(i_shank,1) && spikes.waveform.pos(i_template,2)<=LFP.layers(i_shank,2)
%         FS.layers(i_template) = 3;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(i_shank,2) && spikes.waveform.pos(i_template,2)<=LFP.layers(i_shank,3)
%         FS.layers(i_template) = 4;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(i_shank,3) && spikes.waveform.pos(i_template,2)<=LFP.layers(i_shank,4)
%         FS.layers(i_template) = 5;
%     elseif spikes.waveform.pos(i_template,2)>LFP.layers(i_shank,4)
%         FS.layers(i_template) = 6;
%     end
% 
% 
% 
% end
% %




FS.baseline = baseline;
FS.binwidth = binwidth;
FS.edges = edges;
FS.Offs = Offs;

% get the spikes stats
%PP.SPK.rpv = spikes.isi_vialation; 
FS.SPK.rpv = spikes.isi_vilation;
FS.SPK.stats = spikes.waveform.stats;
FS.SPK.waveform = spikes.waveform.big; 
% also get the spike position
FS.SPK.pos = spikes.waveform.pos; 
%% plot the check new variables
idx = FS.stats_oriN(:,1)==1; 
a=mean(FS.response.first.area_mean(idx,1));
ae=std(FS.response.first.area_mean(idx,1));
b=mean(FS.response.adpt.area_mean(idx,1));
be = std(FS.response.adpt.area_mean(idx,1));

errorbar([1 2],[a b],[ae be])
xlim([0 3])
%%
save('FSN.mat','FS')
%% get the glm fit with increasing number of cells, only for dataset where more than 16 cells
% try all cells first  
% from 8,11,16,24,

rng(0)
nc_sample = [8,11,16,24,30];
nboot = 100;
resp_mat = FS.glm.resp;  % trials x cells (:,FS.stats_Ffit)

targetTrInd = FS.glm.TrlOut'; %vector of zeros and ones that correspond to distractor/target in resp_zs

nc = size(resp_mat,2);
FS.glm.all.pctCorrect = nan(length(nc_sample),nboot);
FS.glm.all.FA = nan(length(nc_sample),nboot);
FS.glm.all.hit = nan(length(nc_sample),nboot);

for iboot = 1:nboot
    for isamp = 1:length(nc_sample)
        if nc < nc_sample(isamp)
            continue
        end
        cell_ind = [];
        cell_ind = randsample(nc,nc_sample(isamp));
        resp_samp = [];
        resp_samp = resp_mat(:,cell_ind);
        %FS.fit.first.prefer(cell_ind)
        [~,~,targetGLM] = glmfit(resp_samp,targetTrInd,'binomial');
        FS.glm.all.weight{isamp}(iboot,:) = targetGLM.beta(2:end);
        %dv_target = mean(targetTrInd);
          dv_target =0.5;  % bias to be very conservative
       [FS.glm.all.pctCorrect(isamp,iboot),FS.glm.all.FA(isamp,iboot),FS.glm.all.hit(isamp,iboot)] = getPctCorr_hoData(...
            resp_samp,targetTrInd,dv_target);
        
    end
end

% store the mean and confidence intervals
FS.glm.all.pctCorrect_mean = nan(length(nc_sample),1);
FS.glm.all.pctCorrect_CI = nan(length(nc_sample),2);
FS.glm.all.FA_mean = nan(length(nc_sample),1);
FS.glm.all.FA_CI = nan(length(nc_sample),2);
FS.glm.all.hit_mean = nan(length(nc_sample),1);
FS.glm.all.hit_CI = nan(length(nc_sample),2);
for isamp = 1:length(nc_sample)
    if nc < nc_sample(isamp)
        continue
    end
    temp = [];
    temp = sort(FS.glm.all.pctCorrect(isamp,:));
    FS.glm.all.pctCorrect_mean(isamp,1) = mean(temp); 
    T_5 = temp(round(0.025*nboot));
    T_95 = temp(round(0.975*nboot));
    FS.glm.all.pctCorrect_CI(isamp,:) = [T_5 T_95];
    temp = [];
    temp = sort(FS.glm.all.FA(isamp,:));
    FS.glm.all.FA_mean(isamp,1) = mean(temp); 
    T_5 = temp(round(0.025*nboot));
    T_95 = temp(round(0.975*nboot));
    FS.glm.all.FA_CI(isamp,:) = [T_5 T_95];
    
    temp = [];
    temp = sort(FS.glm.all.hit(isamp,:));
    FS.glm.all.hit_mean(isamp,1) = mean(temp); 
    T_5 = temp(round(0.025*nboot));
    T_95 = temp(round(0.975*nboot));
    FS.glm.all.hit_CI(isamp,:) = [T_5 T_95];
    
    
end
FS.glm.nc_sample = nc_sample;
%% plot to check 
figure
supertitle('all cells')
subplot(2,2,1)
errorbar(FS.glm.nc_sample,FS.glm.all.FA_mean,abs(FS.glm.all.FA_mean-FS.glm.all.FA_CI(:,1)),abs(FS.glm.all.FA_CI(:,2)-FS.glm.all.FA_mean),'k','Marker','o','LineStyle','none')

xlabel('cell#')
ylabel('FA rate')
xlim([0 40])
ylim([0 1])
axis square 
subplot(2,2,2)
errorbar(FS.glm.nc_sample,FS.glm.all.hit_mean,abs(FS.glm.all.hit_mean-FS.glm.all.hit_CI(:,1)),abs(FS.glm.all.hit_CI(:,2)-FS.glm.all.hit_mean),'k','Marker','o','LineStyle','none')

xlabel('cell#')
ylabel('Hit (90deg)')
xlim([0 40])
ylim([0 1])
axis square 

subplot(2,2,3)
errorbar(FS.glm.nc_sample,FS.glm.all.pctCorrect_mean,abs(FS.glm.all.pctCorrect_mean-FS.glm.all.pctCorrect_CI(:,1)),abs(FS.glm.all.pctCorrect_CI(:,2)-FS.glm.all.pctCorrect_mean),'k','Marker','o','LineStyle','none')

xlabel('cell#')
ylabel('P_correct')
xlim([0 40])
ylim([0 1])
axis square 

%% fit only with well fit neurons and responsive neurons
rng(0)
nc_sample = [8,11,16,24,30];
nboot = 100;
resp_mat = FS.glm.resp(:,FS.stats_F&FS.stats_Ffit);  % trials x cells 

targetTrInd = FS.glm.TrlOut'; %vector of zeros and ones that correspond to distractor/target in resp_zs

nc = size(resp_mat,2);
FS.glm.fit.pctCorrect = nan(length(nc_sample),nboot);
FS.glm.fit.FA = nan(length(nc_sample),nboot);
FS.glm.fit.hit = nan(length(nc_sample),nboot);

for iboot = 1:nboot
    for isamp = 1:length(nc_sample)
        if nc < nc_sample(isamp)
            continue
        end
        cell_ind = [];
        cell_ind = randsample(nc,nc_sample(isamp));
        resp_samp = [];
        resp_samp = resp_mat(:,cell_ind);
        %FS.fit.first.prefer(cell_ind)
%         [~,~,targetGLM] = glmfit(resp_samp,targetTrInd,'binomial');
%         targetWeight = targetGLM.beta(2:end);
%         dv_target = mean(targetTrInd);
          dv_target =0.5;  % bias to be very conservative
       [FS.glm.fit.pctCorrect(isamp,iboot),FS.glm.fit.FA(isamp,iboot),FS.glm.fit.hit(isamp,iboot)] = getPctCorr_hoData(...
            resp_samp,targetTrInd,dv_target);
        
    end
end

% store the mean and confidence intervals
FS.glm.fit.pctCorrect_mean = nan(length(nc_sample),1);
FS.glm.fit.pctCorrect_CI = nan(length(nc_sample),2);
FS.glm.fit.FA_mean = nan(length(nc_sample),1);
FS.glm.fit.FA_CI = nan(length(nc_sample),2);
FS.glm.fit.hit_mean = nan(length(nc_sample),1);
FS.glm.fit.hit_CI = nan(length(nc_sample),2);
for isamp = 1:length(nc_sample)
    if nc < nc_sample(isamp)
        continue
    end
    temp = [];
    temp = sort(FS.glm.fit.pctCorrect(isamp,:));
    FS.glm.fit.pctCorrect_mean(isamp,1) = mean(temp); 
    T_5 = temp(round(0.025*nboot));
    T_95 = temp(round(0.975*nboot));
    FS.glm.fit.pctCorrect_CI(isamp,:) = [T_5 T_95];
    temp = [];
    temp = sort(FS.glm.fit.FA(isamp,:));
    FS.glm.fit.FA_mean(isamp,1) = mean(temp); 
    T_5 = temp(round(0.025*nboot));
    T_95 = temp(round(0.975*nboot));
    FS.glm.fit.FA_CI(isamp,:) = [T_5 T_95];
    
    temp = [];
    temp = sort(FS.glm.fit.hit(isamp,:));
    FS.glm.fit.hit_mean(isamp,1) = mean(temp); 
    T_5 = temp(round(0.025*nboot));
    T_95 = temp(round(0.975*nboot));
    FS.glm.fit.hit_CI(isamp,:) = [T_5 T_95];
    
    
end
%%
figure
supertitle('Good fit')
subplot(2,2,1)
errorbar(FS.glm.nc_sample,FS.glm.fit.FA_mean,abs(FS.glm.fit.FA_mean-FS.glm.fit.FA_CI(:,1)),abs(FS.glm.fit.FA_CI(:,2)-FS.glm.fit.FA_mean),'k','Marker','o','LineStyle','none')

xlabel('cell#')
ylabel('FA rate')
xlim([0 40])
ylim([0 1])
axis square 
subplot(2,2,2)
errorbar(FS.glm.nc_sample,FS.glm.fit.hit_mean,abs(FS.glm.fit.hit_mean-FS.glm.fit.hit_CI(:,1)),abs(FS.glm.fit.hit_CI(:,2)-FS.glm.fit.hit_mean),'k','Marker','o','LineStyle','none')

xlabel('cell#')
ylabel('Hit (90deg)')
xlim([0 40])
ylim([0 1])
axis square 

subplot(2,2,3)
errorbar(FS.glm.nc_sample,FS.glm.fit.pctCorrect_mean,abs(FS.glm.fit.pctCorrect_mean-FS.glm.fit.pctCorrect_CI(:,1)),abs(FS.glm.fit.pctCorrect_CI(:,2)-FS.glm.fit.pctCorrect_mean),'k','Marker','o','LineStyle','none')

xlabel('cell#')
ylabel('P_correct')
xlim([0 40])
ylim([0 1])
axis square 

%%

end 
