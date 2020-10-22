function [spon]=spon(input,Trigger,spikes)

%%
Fs=double(Trigger.SamplingFreq);

Led_N  =  Trigger.data(Trigger.LedID,:)./max(Trigger.data(Trigger.LedID,:));
Led_threshold=Led_N >=0.3;
idx_LED = find(Led_threshold ==1);
diff_LED = diff(idx_LED);

aa=find(diff_LED>=(input.itiTimeMs./1000)*Fs)+1;

LED_time = NaN(length(aa)+1,2);
LED_time(1,1) = idx_LED(1);
LED_time(2:end,1) = idx_LED(aa);
LED_time(1:(end-1),2) = idx_LED(aa-1);
LED_time(end,2) = idx_LED(end);
LED_duration =  mean(diff(LED_time,1,2)./Fs); % in seconds

% add a new baseline value for 2 s before the labjack trigger

Labjack_N  =  Trigger.data(Trigger.labjackID,:)./max(Trigger.data(Trigger.labjackID,:));
Lab_thresh = Labjack_N>=0.3; 
idx_lab = find(Lab_thresh==1);
diff_lab = diff(idx_lab);
bb = find(diff_lab>=(input.itiTimeMs./1000)*Fs)+1;
start_time = NaN(length(bb)+1,1);
start_time(1,1) = idx_lab(1);
start_time(2:end,1) = idx_lab(bb);
% add another baseline value for 2s before the trials immdediately after
% led stimulation
block2 = double(cell2mat(input.tBlock2TrialNumber));
idx_temp = find(diff(block2)==-1)+1;
post_led = start_time(idx_temp);

baseline=2;
binwidth = 0.1;
edges = -baseline:binwidth:6;
basewindow = edges<0;
LED_window = edges>=0& edges<=LED_duration;
edges_off = -5:binwidth:4;% iti:4s 


for i_template = 1:length(spikes.chID) % read the clusters from the spyking circus
    SpikeTimeS = [];
    SpikeTimeS = spikes.spiketime_S{i_template,1}; % already in seconds...
    
    for i_trial = 1:size(LED_time,1)
        
        VstartTimes =[];
        VstartTimes = LED_time(i_trial,1)./Fs;
        spike =[];
        spike_aligned = [];
       
        spike = SpikeTimeS - VstartTimes;
        % arbitually chose duration of 350ms
        spike_aligned = spike(spike>=-baseline & spike<=6); % get a full length that cover longest intervals
        if ~isempty(spike_aligned)
        spike_aligned_new = reshape(spike_aligned,[],length(spike_aligned));      
        else
        spike_aligned_new = NaN;   
        end
        Data{i_template,i_trial} = spike_aligned_new;
        counts_per_bin=[];
        
        counts_per_bin = histc(spike_aligned_new,edges);
        psth.raw(i_template,i_trial,:) = counts_per_bin;
        psth.basebin(i_template,i_trial) = mean(counts_per_bin(basewindow));% the baseline bin used to normalized
        psth.ledbin(i_template,i_trial) = mean(counts_per_bin(LED_window)); % the mean response during LED stim
        
    end
    
    for i_trial = 1:size(start_time,1)
        VstartTimes =[];
        VstartTimes = start_time(i_trial,1)./Fs;
        spike =[];
        spike_aligned = [];
       
        spike = SpikeTimeS - VstartTimes;
        % arbitually chose duration of 350ms
        spike_aligned = spike(spike>=-baseline & spike<=6); % get a full length that cover longest intervals
        if ~isempty(spike_aligned)
        spike_aligned_new = reshape(spike_aligned,[],length(spike_aligned));      
        else
        spike_aligned_new = NaN;   
        end
        counts_per_bin=[];        
        counts_per_bin = histc(spike_aligned_new,edges);
        psth.all.basebin(i_template,i_trial) = mean(counts_per_bin(basewindow));
        
    end 
    
    for i_trial = 1:size(post_led,1)
        VstartTimes =[];
        VstartTimes = post_led(i_trial,1)./Fs;
        spike =[];
        spike_aligned = [];
       
        spike = SpikeTimeS - VstartTimes;
        % arbitually chose duration of 350ms
        spike_aligned = spike(spike>=-baseline & spike<=6); % get a full length that cover longest intervals
        if ~isempty(spike_aligned)
        spike_aligned_new = reshape(spike_aligned,[],length(spike_aligned));      
        else
        spike_aligned_new = NaN;   
        end
        counts_per_bin=[];        
        counts_per_bin = histc(spike_aligned_new,edges);
        psth.postled.basebin(i_template,i_trial) = mean(counts_per_bin(basewindow));
    end 
    % statistic test of whether it is excited or inhibited by led
    % stimulation
    stattest.excite(i_template) = ttest(psth.ledbin(i_template,:),psth.basebin(i_template,:),'Tail','right');
    stattest.inhibit(i_template) = ttest(psth.ledbin(i_template,:),psth.basebin(i_template,:),'Tail','left');
    stattest.pchange(i_template) = (mean(psth.ledbin(i_template,:))-mean(psth.basebin(i_template,:)))./mean(psth.basebin(i_template,:));
end
% align to the led offset to check led recover
for i_template = 1:length(spikes.chID) % read the clusters from the spyking circus
    SpikeTimeS = [];
    SpikeTimeS = spikes.spiketime_S{i_template,1}; % already in seconds...
    
    for i_trial = 1:size(LED_time,1)
        
        VstartTimes =[];
        VstartTimes = LED_time(i_trial,2)./Fs;% align to LED off
        spike =[];
        spike_aligned = [];
        spike = SpikeTimeS - VstartTimes;
        % arbitually chose duration of 350ms
        spike_aligned = spike(spike>=-5 & spike<=4); % get a full length that cover longest intervals
        if ~isempty(spike_aligned)
        spike_aligned_new = reshape(spike_aligned,[],length(spike_aligned));
        end
        Data_off{i_template,i_trial} = spike_aligned_new;
        counts_per_bin=[];
        
        counts_per_bin = histc(spike_aligned_new,edges_off);
        psth.off.raw(i_template,i_trial,:) = counts_per_bin;
        
        
    end
    
end




%plot to see

 PST = [-2 6];
PST = [PST(1)-.0001 PST(2)+.0001];
singleshank=0;
if singleshank==1

 shankchannel=[2,9,3,8,4,7,14,6,5,13,10,1,11,12,15,16 ; 31 19 29 20 28 26 27 24 25 23 30 32 22 21 18 17];
 shanks = [1,2];
 electrodes=16;
else
  shankchannel=[7 14 4 6 8 5 3 13;11 1 12 10 15 2 16 9;32 22 30 21 31 18 19 17;27 26 24 28 25 20 23 29];
  shanks = [1:4];
  electrodes=8;
end
position = reshape(1:32,4,8)';

 figure
for shank=shanks
   
    for j=1:electrodes;
    elec=shankchannel(shank,j);
    idx = find(spikes.chID==elec);
    if ~isempty(idx)
     for i = 1:length(idx)
     subplot(electrodes,length(shanks),position(j,shank))   
     title(['Ch: ' num2str(elec)])
     plotSpikeRaster(Data(idx(i),:),'PlotType','vertline','XLimForCell',PST);
     hold on
     end 
     xlim([-2 6]) 
     hold on 
     vline(0,'b')
     vline( mean(diff(LED_time,1,2)./Fs),'b')
    end
    end
   
end
figure
psth.mean =squeeze( mean(psth.raw,2));
psth.sem = squeeze(std(psth.raw,[],2))./sqrt(size(psth.raw,2));
psth.off.mean =squeeze( mean(psth.off.raw,2));
psth.off.sem = squeeze(std(psth.off.raw,[],2))./sqrt(size(psth.off.raw,2));


for shank=shanks
   
    for j=1:electrodes;
    elec=shankchannel(shank,j);
    idx = find(spikes.chID==elec);
    
    if ~isempty(idx)
     for i = 1:length(idx)
     subplot(electrodes,length(shanks),position(j,shank))   
     
     shadedErrorBar_ch(edges,psth.mean(idx(i),:),psth.sem(idx(i),:));
     title(['Ch: ' num2str(elec)])
     hold on
     end 
     xlim([-2 6])
     
     hold on 
     vline(0,'b')
     vline( mean(diff(LED_time,1,2)./Fs),'b')
    end
    end
   
end

% save the important variables 
spon.edges = edges;
spon.offedges = edges_off;
spon.LED_duration = LED_duration;
spon.binwidth = binwidth;
spon.baseline = baseline;
spon.Data = Data;
spon.spikes = spikes;
spon.Data_off = Data_off;
spon.psth = psth;
spon.stattest = stattest;





end