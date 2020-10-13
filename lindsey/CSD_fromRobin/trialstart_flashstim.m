function [trialstart_s] = trialstart_flashstim(Trigger)

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

% clear up the photodiode trigger and sort the trial start
index_lab=find(labjack_threshold==1);
photo_clean=zeros(1,length(photo_threshold));
photo_clean(index_lab(1):(index_lab(end))) = photo_threshold(index_lab(1):(index_lab(end)));

index_photo=find(photo_clean==1);
diff_index=diff(index_photo);
a=find(diff_index>=(4000./1000)*Fs)+1;

TrialNum=numel(a)+1;

trialstart=zeros(1,TrialNum);
trialstart(1)=index_photo(1);
trialstart(2:end)=index_photo(a(1:(TrialNum-1))); % trial start timestamp align with photodiode signal
trialstart_s = trialstart./Fs;
end