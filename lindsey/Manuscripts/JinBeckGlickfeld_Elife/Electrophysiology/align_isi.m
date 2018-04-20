

function   [Data,psth,Data_full,psth_full]=align_isi(spikes,baseline,binwidth,Time_start_new)
%remember to compsensate according to cycle... store as Data{i_led,1}{i_deg,1}{i_electrode,1}{i_cycle,1}
global Fs 
edges = -baseline:binwidth:0.35;
basebins = round(baseline./binwidth);
Ix_window = edges>=0 & edges <= 0.12; % for the timewindow 120ms after visual start
on_edges =  edges(Ix_window);
edges_full = -baseline:binwidth:5;


%% get the raster data into trial to trial basis

for i_template = 1:length(spikes.chID) % read the clusters from the spyking circus
    SpikeTimeS = [];
    SpikeTimeS = spikes.spiketime_S{i_template,1}; % already in seconds...
    
    for i_trial = 1:length(Time_start_new)
        if length(Time_start_new{i_trial,1})==2 % if there is pair pulse condition
            VstartTimes =[];
            VstartTimes = Time_start_new{i_trial,1}(1)./Fs;
            spike =[];
            spike_aligned = [];
            spike = SpikeTimeS - VstartTimes;
            % arbitually chose duration of 350ms
            spike_aligned = spike(spike>=-baseline & spike<=5); % get a full length that cover longest intervals
            
            Data_full{i_template,1}{i_trial,1} = spike_aligned;
            counts_per_bin=[];
            
            counts_per_bin = histc(spike_aligned,edges_full);
            psth_full.raw{i_template,1}(i_trial,:) = counts_per_bin;
        else
            Data_full = [];
            psth_full = [];
      
        end
        
        for i_cycle=1:length(Time_start_new{i_trial,1})
            
            VstartTimes =[];
            VstartTimes = Time_start_new{i_trial,1}(i_cycle)./Fs;
            spike =[];
            spike_aligned = [];
            spike = SpikeTimeS - VstartTimes;
            % arbitually chose duration of 350ms
            spike_aligned = spike(spike>=-baseline & spike<=0.35);
            
            Data{i_template,1}{i_trial,1}{i_cycle,1} = spike_aligned;
            counts_per_bin=[];
            
            counts_per_bin = histc(spike_aligned,edges);
            psth.raw{i_template,1}{i_trial,1}(i_cycle,:) = counts_per_bin; %
            
            M =[];
            I = [];
            [M,I] = max(counts_per_bin(Ix_window));
            
            if  M~=0
                psth.peak{i_template,1}{i_trial,1}(i_cycle)=M;
                psth.peak_latency_MS{i_template,1}{i_trial,1}(i_cycle) = on_edges(I)*1000;
            else
                psth.peak{i_template,1}{i_trial,1}(i_cycle)=NaN;
                psth.peak_latency_MS{i_template,1}{i_trial,1}(i_cycle) = NaN;
            end
            
            psth.area{i_template,1}{i_trial,1}(i_cycle) = mean(counts_per_bin(Ix_window)); %
            % calculate the BS FR and winodw for ttest analysis
            if i_cycle==1
                psth.baseline{i_template,1}(i_trial,1) = mean(counts_per_bin(edges>=-0.12 & edges <=0)); % mean BS per bin
                psth.First{i_template,1}(i_trial,1) = mean(counts_per_bin(edges>=0 & edges <= 0.12));
            end
            
            if i_cycle==length(Time_start_new{i_trial,1})
                psth.target_baseline{i_template,1}(i_trial,1) = mean(counts_per_bin(edges>=-0.12 & edges <=0));
                psth.target{i_template,1}(i_trial,1)= mean(counts_per_bin(edges>=0 & edges <= 0.12));
            end
            
        end
    end
end












end