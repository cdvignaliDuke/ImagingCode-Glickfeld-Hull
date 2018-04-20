function [LFP] = LFP_layer(NS4,Trigger,input) % or use raw NS6 data
%%
Fs=double(Trigger.SamplingFreq);
 %[Time_start_new,Time_stop_new,photo_duration_MS,photo_inter_MS,~]=PhotoISI_nolabjack(Trigger,input);
[Time_start_new,Time_stop_new,photo_duration_MS,photo_inter_MS,~]=PhotoISI(Trigger,input);
%%
baseline=0.2;
if isfield(input,'tStimOffTimes')
    StimOffs = double(cell2mat(input.tStimOffTimes));
else
    StimOffs = double(cell2mat(input.tStimOffTimeMs));
    
end

if length(unique(StimOffs))>3  % for pair pulse stimulus only
LFP_SF = NS4.MetaTags.SamplingFreq;
baselinebins = baseline*LFP_SF;
durationbins = (0.1+0.4)*LFP_SF; % get extra 400ms after visual onset

LFP_Stamp = NS4.MetaTags.Timestamp
edges = [-baseline:1/LFP_SF:(0.1+0.4)];

Data_filter = {};
Data_combine = {};

% only get the LFP data for first pulse with ISI>250 and last pulse isi4s


for Elec_Num = 1:length(NS4.MetaTags.ChannelID) % read the clusters from the spyking circus
    Elec_num = NS4.MetaTags.ChannelID(Elec_Num);
    LFP_temp = [];
    LFP_temp = double(NS4.Data(Elec_Num,:));
    LFP_Filter = [];
    LFP_Filter = LFPfilter (LFP_temp,[0.1 200], LFP_SF); % matches the filter in haider's paper,but cannot be lower htan 1..
    i=1;
    for i_trial = 1:size(Time_start_new,1)
        if StimOffs(i_trial)>250 &&  StimOffs(i_trial)<4000
            VstartTimes =[];
            VstartTimes = Time_start_new{i_trial,1}(1,1)./Fs;
            LFP.data(Elec_num,i,:) = LFP_Filter ( (int64(VstartTimes*LFP_SF)-baselinebins):(int64(VstartTimes*LFP_SF)+durationbins));
            
            i=i+1;
        end
        
        if  StimOffs(i_trial)==4000
            VstartTimes =[];
            VstartTimes = Time_start_new{i_trial,1}(1,1)./Fs;
            LFP.data(Elec_num,i,:) = LFP_Filter ( (int64(VstartTimes*LFP_SF)-baselinebins):(int64(VstartTimes*LFP_SF)+durationbins));
            
            i=i+1;
            VstartTimes =[];
            VstartTimes = Time_start_new{i_trial,1}(1,2)./Fs;
            LFP.data(Elec_num,i,:) = LFP_Filter ( (int64(VstartTimes*LFP_SF)-baselinebins):(int64(VstartTimes*LFP_SF)+durationbins));
            
            i=i+1;
            
        end
    end
    
    
    
    
end

elseif length(unique(StimOffs))==3% this is for randomISI, only get the response to the first stimulus 
    
LFP_SF = NS4.MetaTags.SamplingFreq;
baselinebins = baseline*LFP_SF;
durationbins = (0.1+0.4)*LFP_SF; % get extra 400ms after visual onset

LFP_Stamp = NS4.MetaTags.Timestamp
edges = [-baseline:1/LFP_SF:(0.1+0.4)];

Data_filter = {};
Data_combine = {};

% only get the LFP data for first pulse with ISI>250 


for Elec_Num = 1:length(NS4.MetaTags.ChannelID) % read the clusters from the spyking circus
    Elec_num = NS4.MetaTags.ChannelID(Elec_Num);
    LFP_temp = [];
    LFP_temp = double(NS4.Data(Elec_Num,:));
    LFP_Filter = [];
    LFP_Filter = LFPfilter (LFP_temp,[0.1 200], LFP_SF); % matches the filter in haider's paper,but cannot be lower htan 1..
   
    i=1;
    for i_trial = 1:size(Time_start_new,1)
        if input.tStimOffTimes{i_trial}(1)>250 
            VstartTimes =[];
            VstartTimes = Time_start_new{i_trial,1}(1,1)./Fs;
            LFP.data(Elec_num,i,:) = LFP_Filter ( (int64(VstartTimes*LFP_SF)-baselinebins):(int64(VstartTimes*LFP_SF)+durationbins));
            
            i=i+1;
        end
        
      
    end
    
    
    
    
end    
    
else  % this is for FS orien, get the response for each cycle
LFP_SF = NS4.MetaTags.SamplingFreq;
baselinebins = baseline*LFP_SF;
ITIbins = LFP_SF; % get one second long ITI responses
durationbins = (0.1+0.4)*LFP_SF; % get extra 400ms after visual onset

edges = [-baseline:1/LFP_SF:(0.1+0.4)];
ITIedges = [0:1/LFP_SF:1]; % baseline responses 1 s before vistim is on
Data_filter = {};
Data_combine = {};

% only get the LFP data for each pulses
ChannelID = [];
ChannelID = NS4.MetaTags.ChannelID(NS4.MetaTags.ChannelID<33);

for Elec_Num = 1:length(ChannelID) % read the clusters from the spyking circus
    Elec_num = ChannelID(Elec_Num);
    LFP_temp = [];
    LFP_temp = double(NS4.Data(Elec_Num,:));
    LFP_Filter = [];
    LFP_Filter = LFPfilter (LFP_temp,[0.1 200], LFP_SF); % matches the filter in haider's paper,but cannot be lower htan 1..[1 200]
   
    for i_trial = 1:size(Time_start_new,1)
        for i_cycle = 1:(input.maxCyclesOn + 1)% it is usually has 6 cycles
            VstartTimes =[];
            VstartTimes = Time_start_new{i_trial,1}(i_cycle)./Fs;
            temp_trl = [];
            temp_trl =  LFP_Filter ( (int64(VstartTimes*LFP_SF)-baselinebins):(int64(VstartTimes*LFP_SF)+durationbins));
            
            LFP.Data(Elec_num,i_trial,i_cycle,:) = temp_trl;
            % zero the baseline for each trial, 50ms epoke before visual
            % onset
            offset = 0;
            offset = mean(temp_trl(1501:2000));
            
            LFP.Data_zerobase(Elec_num,i_trial,i_cycle,:) = temp_trl - offset; 
            
            
            % get the first cycle responses for calculating layers
            if i_cycle ==1
                LFP.data(Elec_num,i_trial,:) = LFP_Filter ( (int64(VstartTimes*LFP_SF)-baselinebins):(int64(VstartTimes*LFP_SF)+durationbins));
                % get the baseline LFP during the ITI
                LFP.base(Elec_num,i_trial,:) = LFP_Filter ( (int64(VstartTimes*LFP_SF)-ITIbins):int64(VstartTimes*LFP_SF));
            end 
            
        end
        
        
    end
    
    
    
    
end



end

 


%% check bad channels and replace it by mean of channels in between, this is for old conditions?


poly=0;
if poly==1
    shankchannel=[2,9,3,8,4,7,14,6,5,13,10,1,11,12,15,16 ; 31 19 29 20 28 26 27 24 25 23 30 32 22 21 18 17];
    e_distance(1,:) = [25:50:775]; % shank1
    e_distance(2,:) = [0:50:750]; % shank2
    
    for shank = 1:size(shankchannel,1)
        figure
        for j = 1:size(shankchannel,2)
            elec=shankchannel(shank,j);
            LFP.mean(shank,j,:) = mean(squeeze(LFP.data(elec,:,:)),1)-mean(mean(squeeze(LFP.data(elec,:,1:baselinebins)),1));
            %         LFP.trial(shank,j,:,:) = squeeze(LFP.data(elec,:,:)) - mean(mean(squeeze(LFP.data(elec,:,:)),1));
            subplot(4,4,j)
            plot(edges,squeeze(LFP.mean(shank,j,:)))
            title(['shank: ', num2str(shank), ' ', 'channel: ', num2str(elec)])
            %ylim([-1500 1500]) % for V1
            ylim([-500 500]) % for HVAs
            xlim([-0.2 0.35])
        end
    end
    
else
    shankchannel=[7 14 4 6 8 5 3 13;11 1 12 10 15 2 16 9;32 22 30 21 31 18 19 17;27 26 24 28 25 20 23 29];
    e_distance = repmat([0:100:700],size(shankchannel,1),1);
    for shank = 1:size(shankchannel,1)
        figure
        for j = 1:size(shankchannel,2)
            elec=shankchannel(shank,j);
            LFP.mean(shank,j,:) = mean(squeeze(LFP.data(elec,:,:)),1)-mean(mean(squeeze(LFP.data(elec,:,1:baselinebins)),1));
            %         LFP.trial(shank,j,:,:) = squeeze(LFP.data(elec,:,:)) - mean(mean(squeeze(LFP.data(elec,:,:)),1));
            subplot(3,3,j)
            plot(edges,squeeze(LFP.mean(shank,j,:)))
            title(['shank: ', num2str(shank), ' ', 'channel: ', num2str(elec)])
            %ylim([-1500 1500]) % for V1
            ylim([-1000 1000]) % for HVAs
        end
    end
    
end

%% store the broken channels
LFP.broken_ch = 22;

%% mannually get rid of broken channels
if poly==1
broken_ch = 22;
shank=2;
j = find(shankchannel(2,:)==22);
% replace it by the mean of surround signal
LFP.mean(shank,j,:) =mean([squeeze(LFP.mean(shank,j+1,:)),squeeze(LFP.mean(shank,j-1,:))],2);
else
broken_ch = 22;
shank=3;
j = find(shankchannel(3,:)==22);
% replace it by the mean of surround signal
LFP.mean(shank,j,:) =mean([squeeze(LFP.mean(shank,j+1,:)),squeeze(LFP.mean(shank,j-1,:))],2);  
    
end 
%% re plot for checking
if poly==1
    for shank = 1:size(shankchannel,1)
        figure
        for j = 1:size(shankchannel,2)
            elec=shankchannel(shank,j);
            
            subplot(4,4,j)
            plot(edges,squeeze(LFP.mean(shank,j,:)))
            title(['shank: ', num2str(shank), ' ', 'channel: ', num2str(elec)])
                     ylim([-800 800])
            %ylim([-500 500]) % for HVAs
             xlim([-0.2 0.35])
        end
    end
    
else
    for shank = 1:size(shankchannel,1)
        figure
        for j = 1:size(shankchannel,2)
            elec=shankchannel(shank,j);
            
            subplot(3,3,j)
            plot(edges,squeeze(LFP.mean(shank,j,:)))
            title(['shank: ', num2str(shank), ' ', 'channel: ', num2str(elec)])
                     ylim([-1500 1500])
           % ylim([-1000 1000]) % for HVAs
        end
    end
    
    
end

% not get rid of noisy trials for now 

%% calculate CSD
% shankchannel=[2,9,3,8,4,7,14,6,5,13,10,1,11,12,15,16 ; 31 19 29 20 28 26 27 24 25 23 30 32 22 21 18 17];
% e_distance(1,:) = [25:50:775]; % shank1
% e_distance(2,:) = [0:50:750]; % shank2
% baseline=0.2;
% edges = [-baseline:1/10000:(0.1+0.4)];

for shank = 1:size(shankchannel,1)
    d=-diff(squeeze(LFP.mean(shank,:,:)),2,1);% approximate second spatial derivative
    d_i = interp2(d); % doubles the samples in every direction
    d_f = d_i(:,1:2:size(d_i,2)); % downsamples the times again
    test = linspace(min(e_distance(shank,:)),max(e_distance(shank,:)),size(d_i,1));
    figure
    pcolor(edges,(-test),d_f);
    xlim([-0.1 0.3])
    shading interp;
    colormap(jet)
    caxis manual
    
     caxis([-100 100]) % for V1 /SC
    %caxis([-100 100]) % for HVAs
    colorbar
    title(['shank:', num2str(shank)])
%     hold on
%     scatter(0.05,-400,'MarkerEdgeColor',[1 1 1])

  

    
end
%% if the code is FSorien, then run the functions for sumarize LFP 
% need to get broken channel fixed
poly = 1;
LFP.sum = LFP_sum(LFP,input, edges, baselinebins,poly,LFP.broken_ch); 

%% get the handle of layers on one of the csd figures
if poly==1
shank=2;
d=-diff(squeeze(LFP.mean(shank,:,:)),2,1);% approximate second spatial derivative
d_i = interp2(d); % doubles the samples in every direction
d_f = d_i(:,1:2:size(d_i,2)); % downsamples the times again
test = linspace(min(e_distance(shank,:)),max(e_distance(shank,:)),size(d_i,1));
figure
pcolor(edges,(-test),d_f);
xlim([-0.1 0.3])
shading interp;
colormap(jet)
caxis manual


 caxis([-200 200]) % for V1 /SC
%caxis([-100 100]) % for HVAs
colorbar
title(['shank:', num2str(shank)])
[x,y] = ginput(4); % for V1 csd, draw line on the boundaries of layers
%[x,y]= ginput(3);% for HVA csd, draw line on the boundaries around two sinks.

% draw the boundary lines
hold on
hline(y,'k:')
vline([0 0.1],'K:')
LFP.shank = shank;
else
% for four shank probe
for shank = 1:4
d=-diff(squeeze(LFP.mean(shank,:,:)),2,1);% approximate second spatial derivative
d_i = interp2(d); % doubles the samples in every direction
d_f = d_i(:,1:2:size(d_i,2)); % downsamples the times again
test = linspace(min(e_distance(shank,:)),max(e_distance(shank,:)),size(d_i,1));
figure
pcolor(edges,(-test),d_f);
xlim([-0.1 0.3])
shading interp;
colormap(jet)
caxis manual


%caxis([-200 200]) % for V1 /SC
caxis([-100 100]) % for HVAs
colorbar


title(['shank:', num2str(shank)])
[x(shank,1:4),y(shank,1:4)] = ginput(4); % for V1 csd, draw line on the boundaries of layers
%[x(shank,1:3),y(shank,1:3)]= ginput(3);
hold on
hline(y(shank,:),'k:')
vline([0 0.1],'K:')
LFP.shank = shank; % which shank is chosen     
    
end


end 

%% save the layers and also LFP
LFP.Data=[];
LFP.data = [];
LFP.Data_zerobase = [];
LFP.base = [];
LFP.mean = [];
LFP.layers = abs(y);
LFP.edges = edges;
LFP.shankchannel = shankchannel;
LFP.e_distance = e_distance;

save('layers.mat','LFP')
%save('LFP.mat','LFP')
%% save the data
LFP.sum = [];
LFP.ITIedges = ITIedges;
% if all layers are 0, suggests that cannot distinguish between layers
save('layers.mat','LFP')

% LFP.data = [];
% LFP.Data = [];

% layer1 when > =layers(1)
% layer2/3 when<layers(1) && when >=layers(2)
% layer4 when <layers(2) && when >=layers(3)
% layer5 when <layers(3) && when>=layers(4)
% layer6 when <layers(4)

% if layer values equals=0, suggests no layers there


end


