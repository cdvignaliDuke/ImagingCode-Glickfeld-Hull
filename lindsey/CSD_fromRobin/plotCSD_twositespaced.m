
%% Load the trigger data and NS4 file 

[trialstart_s] = trialstart_twostim(Trigger); %will need to change based on experiment used (two stim vs. flashing stim)
%[trialstart_s] = trialstart_flashstim(Trigger);

Fs=double(Trigger.SamplingFreq);
baseline=0.2;

LFP_SF = NS4.MetaTags.SamplingFreq;
baselinebins = baseline*LFP_SF;
durationbins = (0.5+1)*LFP_SF; % get extra 400ms after visual onset

LFP_Stamp = NS4.MetaTags.Timestamp
edges = [-baseline:1/LFP_SF:(0.5+1)]; 

Data_filter = {};
Data_combine = {};

%% Concatenate across timestamps if data is divided up (don't know why this is happening) 

if size(NS4.Data,2)>1 
    for i = 2:size(NS4.Data,2) 
    NS4.Data{1,1} = horzcat(NS4.Data{1,1}, NS4.Data{1,i});
    end 
end 

if size(NS4.Data,2)>1 
NS4.Data = NS4.Data{1,1}
end

%% Low pass filter the NS4 file data 

for Elec_Num = 1:length(NS4.MetaTags.ChannelID) % read the clusters from the spyking circus
    Elec_num = NS4.MetaTags.ChannelID(Elec_Num);
    LFP_temp = [];
    LFP_temp = double(NS4.Data(Elec_Num,:));
    LFP_Filter = [];
    LFP_Filter = LFPfilter (LFP_temp,[0.1 200], LFP_SF); % matches the filter in haider's paper,but cannot be lower than 1..
    i=1;
    for i_trial = 1:length(trialstart_s) 
            VstartTimes =[];
            VstartTimes = trialstart_s(i_trial);
            LFP.data(Elec_num,i_trial,:) = LFP_Filter((int64((VstartTimes*LFP_SF)-baselinebins):(int64((VstartTimes*LFP_SF)+durationbins))));
    end
end 

%% Generate mean LFP for each channel 

shankchannel=[2,9,3,8,4,7,14,6,5,13,10,1,11,12,15,16 ; 31 19 29 20 28 26 27 24 25 23 30 32 22 21 18 17];
e_distance(1,:) = [25:50:775]; % shank1
e_distance(2,:) = [0:50:750]; % shank2
for shank = 1:size(shankchannel,1)
    figure
    for j = 1:size(shankchannel,2)
        elec=shankchannel(shank,j);
        LFP.mean(shank,j,:) = mean(squeeze(LFP.data(elec,:,:)),1)-mean(mean(squeeze(LFP.data(elec,:,1:baselinebins)),1));
        subplot(4,4,j)
        plot(edges,squeeze(LFP.mean(shank,j,:)));
        title(['shank: ', num2str(shank), ' ', 'channel: ', num2str(elec)])
        ylim([-1500 500])
        xlim([-0.2 0.35])
    end
end

%% check the channels and replace with mean of surrounding channels if noisy 

LFP_sep_s1 = bsxfun(@plus,(100:100:16*100)',squeeze(LFP.mean(1,[1:16],:))); %adds 100 to each channel 
figure; plot(edges, LFP_sep_s1) 
title('Mean stimulus-evoked LFP (column 1)')
xlabel('Time(s)')
ylabel('Mean LFP (with arbitrary 100uVolt offset)') 

LFP_sep_s2 = bsxfun(@plus,(100:100:16*100)',squeeze(LFP.mean(2,[1:16],:))); %adds 100 to each channel 
figure; plot(edges, LFP_sep_s2) 
title('Mean stimulus-evoked LFP (column 2)')
xlabel('Time(s)')
ylabel('Mean LFP (with arbitrary 100uVolt offset)') 

% These are all commonly noisy channels 
LFP.mean(2,13,:) = mean([squeeze(LFP.mean(2,14,:)), squeeze(LFP.mean(2,12,:))], 2);
LFP.mean(1,10,:) = mean([squeeze(LFP.mean(1,9,:)), squeeze(LFP.mean(1,11,:))], 2);
%LFP.mean(1,11,:) = mean([squeeze(LFP.mean(1,10,:)), squeeze(LFP.mean(1,12,:))], 2);

%% Normalize LFP data to max sink 

for i = 1:size(LFP.mean,2)
    LFP.mean_norm{1,1}(i,:) = squeeze(LFP.mean(1,i,:));
end 

min_LFP = min(LFP.mean_norm{1,1}(:)); 

find(LFP.mean_norm{1,1}(:) == -1);

LFP.mean_norm{1,1} = -(LFP.mean_norm{1,1}./min_LFP);

LFP_sep_s2 = bsxfun(@plus,(.2:.2:16*.2)',squeeze(LFP.mean_norm{1,1}([1:16],:))); %adds 100 to each channel 
figure; plot(edges, LFP_sep_s2) 
title('Mean stimulus-evoked LFP (column 2)')
xlabel('Time(s)')
ylabel('Mean LFP (with arbitrary 100uVolt offset)') 


%% Plot the CSD

e_distance(1,:) = [25:50:775]; % shank1
e_distance(2,:) = [0:50:750]; % shank2
baseline=0.2;
edges = [-baseline:1/10000:(0.5+1)];

d=-diff(squeeze(LFP.mean(1,1:2:16,:)),2,1);% approximate second spatial derivative for every other channel 
min_second_deriv = min(d(:)); % find min
d_norm = -(d./min_second_deriv); %normalize second spatial derivative to mib
d_i = interp2(d_norm); % doubles the samples in every direction
d_f = d_i(:,1:2:size(d_i,2)); % downsamples the times again
test = linspace(e_distance(1,1),e_distance(1,15),size(d_i,1));
figure
pcolor(edges,(-test),d_f);
xlim([-0.1 1.5])
shading interp;
colormap(jet)
caxis manual
caxis([-1 1])
colorbar

