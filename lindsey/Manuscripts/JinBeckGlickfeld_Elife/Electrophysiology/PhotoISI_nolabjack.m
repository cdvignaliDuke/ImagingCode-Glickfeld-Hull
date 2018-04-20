function [Time_start_new,Time_stop_new,photo_duration_MS,photo_inter_MS,LED_time]=PhotoISI_nolabjack(Trigger,input)
global Fs stimOnTimeS Offs


target_compensate = 0;


% define the compensate adjustment for photodiode signal based on elevation
if unique(double(cell2mat(input.tGratingHeightDeg)))==30
    compensate =8-unique(double(cell2mat(input.tGratingElevationDeg)))*0.4;
else
    compensate =0;
end

Cycle = double(unique(cell2mat(input.tCyclesOn))); % number of cycles
Cycle_trl = double(cell2mat(input.tCyclesOn)); % how many cycles each trial
Trialnum = input.trialSinceReset;

stimOnTimeS = double(input.stimOnTimeMs)./1000;

% try to align the onset trigger times trial to trial basis?
if isfield(input,'tStimOffTimes') 
Offs = unique(double(cell2mat(input.tStimOffTimes))); % how many off conditions
else
Offs = unique(double(cell2mat(input.tStimOffTimeMs)));
end




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
photo_threshold = photo_filter_N>=0.2;



% labjack_N  =  Trigger.data(Trigger.labjackID,:)./max(Trigger.data(Trigger.labjackID,:));
% labjack_threshold=labjack_N >=0.5;
% without labjack trigger, get the start and end trials
e=find(photo_filter_N>=0.2);
f=diff(e);
g = find(f>(double(input.itiTimeMs)./1000)*Fs); % uesed to get the start 
h = find(f>5000); 


if input.doBlock2 ==1
if  isfield(input, 'tStimOffTimes')   
input.tStimOffTimes = cellfun(@(x) x(2:end),input.tStimOffTimes, 'UniformOutput', false);
end
% also get the led trigger 
Led_N  =  Trigger.data(Trigger.LedID,:)./max(Trigger.data(Trigger.LedID,:));
Led_threshold=Led_N >=0.3;

end

%% clear up the photodiode trigger and sort the trial start
% index_lab=find(labjack_threshold==1);
photo_clean=zeros(1,length(photo_threshold));
% photo_clean(e(g(1)+1)-1000:(e(h(end)+1)-1000))=photo_threshold(e(g(1)+1)-1000:(e(h(end)+1)-1000));
photo_clean(e(g(1)+1)-1000:(e(h(end)+1)-1000))=photo_threshold(e(g(1)+1)-1000:(e(h(end)+1)-1000)); % for PP

index_photo=find(photo_clean==1);
diff_index=diff(index_photo);
a=find(diff_index>=(input.itiTimeMs./1000)*Fs)+1; 

TrialNum=numel(a)+1;
%%------manually put trial number if crashes--------
%  TrialNum = 66;
if TrialNum <Trialnum
   manual=1;
else
    manual=0;
end
TrialNum=min([TrialNum,Trialnum]);

%% manually modify TrialNum if code breaks!
% need to modify the corresponding Index !
% Index_part = {}; % store as Index{Led,1}{Orien,1} orientation level start from smaller to bigger..without led to with led
% for i_led = 1:2
%     for orien_N = 1:numel(Orienlevel)
%         Index_part{i_led,1}{orien_N,1} = intersect(Index{i_led,1}{orien_N,1},1:TrialNum);
%     end
% end

trialstart=zeros(1,TrialNum);
trialstart(1)=index_photo(1);
trialstart(2:end)=index_photo(a(1:(TrialNum-1))); % trial start timestamp align with photodiode signal


%% get the timestamp of each cycle start

if input.doBlock2 ==1
    Led = cell(TrialNum,1);
end

photo_cycle_start = cell(TrialNum,1);
photo_cycle_stop = cell(TrialNum,1);
photo = cell(TrialNum,1);
Time_start = cell(TrialNum,1);
Time_stop = cell(TrialNum,1);
Time_start_new = cell(TrialNum,1);
Time_stop_new = cell(TrialNum,1);
photo_duration = cell(TrialNum,1);
photo_inter = cell(TrialNum,1);
photo_duration_MS = cell(TrialNum,1);
photo_inter_MS =cell(TrialNum,1);
for i=1:TrialNum
    Time_start{i,1}=zeros(1,(Cycle_trl(i)+1));
    Time_stop{i,1}=zeros(1,(Cycle_trl(i)+1));
    Time_start_new{i,1} = zeros(1,(Cycle_trl(i)+1));
    Time_stop_new{i,1} = zeros(1,(Cycle_trl(i)+1));
    photo_cycle_start{i,1}=zeros(1,(Cycle_trl(i)+1)); % the last one always is
    photo_cycle_stop{i,1}=zeros(1,(Cycle_trl(i)+1));
    
    if isfield(input, 'tStimOffTimes')
    offs_trl = double(input.tStimOffTimes{i})./1000; % off times for each trial, the second last one is the one before target on
    duration=int64(Fs.*(sum(offs_trl)+stimOnTimeS*(Cycle_trl(i)+1))); % for the duration of each trial
    else
    offs_trl = double(input.tStimOffTimeMs{i})./1000; % off times for each trial, the second last one is the one before target on
    duration=int64(Fs.*(offs_trl+stimOnTimeS)*(Cycle_trl(i)+1)); % for the duration of each trial     
    end 
   
    
    try
        % 3000 is the baseline...
        photo{i,:}=photo_clean(trialstart(i)-3000:trialstart(i)+duration); %use arryfun later
        
    catch
        warning('trial cycles not complete');
        % 3000 is the baseline...
        photo{i,:}=zeros(1,3000+duration);
    end
    
     if input.doBlock2 ==1
     try
        Led{i,:}=Led_threshold(trialstart(i)-3000:trialstart(i)+duration); % 3000 is the baseline...
        
        
    catch
        warning('trial cycles not complete');
        Led{i,:}=zeros(1,3000+duration); % 3000 is the baseline...
      
    end    
        
     end
    
    index_cycle=[];
    diff_cycle=[];
    c=[];
    index_cycle=find(photo{i,:}==1);
    diff_cycle=diff(index_cycle);
    if length(Offs)>2 % deal with pair pulse condition
    
    c=intersect(find(diff_cycle>=6000),find(diff_cycle<=4*Fs+6000))+1;
        
    else
%     c=intersect(find(diff_cycle>=6000),find(diff_cycle<=52000))+1;

    c=intersect(find(diff_cycle>=6000),find(diff_cycle<=52000))+1;
    end
    
    c_new = c(1:Cycle_trl(i));
    photo_cycle_start{i,1}(1)= index_cycle(1);
    photo_cycle_start{i,1}(2:(Cycle_trl(i)+1))=index_cycle(c_new);
    photo_cycle_stop{i,1}(1:Cycle_trl(i))=index_cycle(c_new-1);
   % modify hard coded in non labjact condition
    %photo_cycle_stop{i,1}((Cycle_trl(i)+1))=index_cycle(c_new)+(index_cycle(c_new-1)-index_cycle(1));
    photo_cycle_stop{i,1}((Cycle_trl(i)+1))=index_cycle(end);
    photo_duration{i,1}=photo_cycle_stop{i,1}-photo_cycle_start{i,1};
    photo_inter{i,1}=photo_cycle_start{i,1}(2:(Cycle_trl(i)+1))-photo_cycle_stop{i,1}(1:Cycle_trl(i));
    photo_duration_MS{i,1}=photo_duration{i,1}.*1000/Fs;
    photo_inter_MS{i,1}=photo_inter{i,1}.*1000./Fs;
    
    trial_duration=0;
    
    for j=1:(Cycle_trl(i)+1)
        Time_start{i,1}(j)=trialstart(i)+trial_duration;
        Time_start_new{i,1}(j) =Time_start{i,1}(j) +(compensate.*Fs./1000); % compensate the difference between photodiode and visual stimus signal
        
        
        
        if j<=Cycle_trl(i)
            trial_duration=trial_duration + photo_duration{i,1}(j)+photo_inter{i,1}(j);
        end
        
        
        Time_stop{i,1}(j)=Time_start{i,1}(j)+photo_duration{i,1}(j);
        Time_stop_new{i,1}(j) = Time_stop{i,1}(j)+(compensate.*Fs./1000);
    end
    
    
    
    
    
    
end
%% get the LED time stamp relative to the visual onset 
if input.doBlock2 ==1
block2 = cell2mat(input.tBlock2TrialNumber); % one is led trials

temp = Led(block2==1,:);
temp1 = photo(block2==1,:);

idx = cellfun(@(x) find(x==1),temp,'UniformOutput',false);
idx1 = cellfun(@(x) find(x==1),temp1,'UniformOutput',false);
time = cell2mat(cellfun(@(x,y) (y(1,1)-x(1,1))./Fs,idx,idx1,'UniformOutput',false));
% calculat led duration
led_dura = cell2mat(cellfun(@(x) sum(x==1)./Fs,temp,'UniformOutput',false));
LED_time(1) = -mean(time);
LED_time(2)= LED_time(1)  + mean(led_dura); % duration is 200ms 
LED_time = LED_time - compensate./1000;
else
    LED_time = NaN;
end







%% check if it is correct sort out non led condition
trial=100;
    
figure

plot(photo{trial,:},'g')
%% check if it is correct

figure
vline(Time_start_new{trial,1},'k')
hold on
vline(Time_stop_new{trial,1},'r')




end