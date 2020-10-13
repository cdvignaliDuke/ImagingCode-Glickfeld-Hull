
function [trialstart_s] = trialstart_twostim(Trigger)
% get the all the normalized and filtered trigger
Fs=double(Trigger.SamplingFreq);
dt = 1/Fs;

fpass  =    [80 400];  % pass freq
fstop  =    [70 500];  % stop freq outside pass
Rpass  =    0.5;       % Attenuation(dB) for passband
Astop  =    30;        % Attenuation(dB) for stopband
n      =    cheb2ord(fpass/Fs*2,fstop/Fs*2,Rpass,Astop);   % order of chey filter

[z,p,k] =   cheby2(n,Astop,fstop/Fs*2);   % zeros, poles, and gain
[s,g]   =   zp2sos(z,p,k);                  % create second order section
Hd      =   dfilt.df2sos(s,g);                 % dfilt object
photo_filter  =    filtfilthd(Hd,Trigger.data(Trigger.photoID,:));    % apply filter


photo_Subtract  = Trigger.data(Trigger.photoID,:)-mean(Trigger.data(Trigger.photoID,:));
photo_N         =   photo_Subtract./max(photo_Subtract);


photo_filter_N  =   photo_filter./max(photo_filter);
photo_threshold = photo_filter_N>=0.5; % threshold photodiode signal from photodiode noise?

%photo_threshold = photo_N  >=0.4;% deal with some odd conditions, date: 161208 ; 161215 ; 161214

labjack_N  =  Trigger.data(Trigger.labjackID,:)./max(Trigger.data(Trigger.labjackID,:));
labjack_threshold=labjack_N >=0.5; % threshold labjack signal from labjack noise? 


% photodiode analysis: align with trigger, plot on an off times 

diff_labjack_thresh = diff(labjack_threshold);
index_lab=find(labjack_threshold==1);
trigger_on = find(diff_labjack_thresh == 1); 
trigger_off = find(diff_labjack_thresh == -1);
trigger_diff = trigger_on(1:length(trigger_off)) - trigger_off(1:length(trigger_off)); 
stimOne_idx = find(trigger_diff < -1000); 
trig_stimOne = trigger_on(stimOne_idx);
photo_clean=zeros(1,length(photo_threshold));
photo_clean(trig_stimOne(1):(trig_stimOne(end)+int64(Fs)))=photo_threshold(trig_stimOne(1):(trig_stimOne(end)+int64(Fs)));

photo_clean(trig_stimOne(1):(trig_stimOne(end)))=photo_threshold(trig_stimOne(1):(trig_stimOne(end)));


photo_on = find(photo_clean == 1);

photo_on_time = 1:length(photo_clean); 

diff_on = diff(photo_on); 
diff_on_gap = find(diff_on > 1000);

trialstart = photo_on(diff_on_gap)

trial_time_total = diff(trialstart)./30000

for i = 1:(length(trialstart)-1)
    indiv_trials{i,1} = photo_clean(trialstart(i):(trialstart(i+1)-1)); 
    indiv_trials{i,1}(2,:) = photo_on_time(trialstart(i):(trialstart(i+1)-1)); 
end 

for i = 1:(length(trialstart)-1)
    photodiode_on_idx = find(indiv_trials{i,1}(1,:) == 1); 
    photodiode_stimon(i) = (indiv_trials{i,1}(2,photodiode_on_idx(2)));
    photodiode_stimon_length(i) = (photodiode_on_idx(end) - photodiode_on_idx(2))./30000;
    photodiode_ITItime(i) = ((indiv_trials{i,1}(2,photodiode_on_idx(2))-1) - (indiv_trials{i,1}(2,photodiode_on_idx(1))+1))./30000; 
end 


index_photo=find(photo_clean==1);
diff_index=diff(index_photo);
a=find(diff_index>=(6000./1000)*Fs)+1; 

TrialNum=numel(a)+1;

trialstart=zeros(1,TrialNum);
trialstart(1)=index_photo(1);
trialstart(2:end)=index_photo(a(1:(TrialNum-1))); % trial start timestamp align with photodiode signal
%index is on sampling rate space so divide by 30000 
trialstart_s = trialstart./30000;

end