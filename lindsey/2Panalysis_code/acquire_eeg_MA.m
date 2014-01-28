%acqures daq from nidaq

clear

%parameters to change
DIR1 = '110630';
MOUSE = '110630_CH1_Y6_CH2_DR5_CH3';
Nchan = 2; %8; # of channels to acquire
Fs = 1000; %sampling; 98304; %512*256*60;
%T=60*60; %sec
T = 120*60; %60; %in seconds, 60*(# of minutes);   15.5 for 8dir x 3 SF; 19.5; %20.5; %19.5; %60*19.5; %; %16.5; %16; %60*1; %20; % %200; %450



drifting = 1;
Ftemp = 2; 
file_add = ['rec_',MOUSE,'inrig_grat_drift',num2str(drifting),'_Fhz',num2str(Ftemp),'_Per_20']
PWD = ['C:\','Documents and Settings','\','Mark Andermann','\Desktop\Andermann_data\'];
PWD_save = [DIR1,'\'];

%PWD =['C:\',''','Documents and Settings',''','\',''','Mark Andermann',''','\Desktop\Andermann_data\Mfiles_daq'];
cd(PWD)
file00 = 'mouse_training';
DIR = dir(PWD_save);
if isempty(DIR)
    eval(['!mkdir ',PWD_save]);
end

a = clock;
a3 = [num2str(a(1)),'-', ...
    num2str(a(2)),'-', ...
   num2str(a(3)),'_', ...
   num2str(a(4)),'-', ...
   num2str(a(5)),'-', ...
   num2str(round(a(6)))];

filesave = [PWD_save,'eegacq_',file00,a3,'_',file_add,'.mat'];
disp(['opening ',filesave])

%set # of channels, etc

T2 = round(T*Fs);

ai=analoginput('nidaq','Dev1');
%set(ai,'timeout',T + 10);
set(ai,'SampleRate',Fs);
set(ai,'Timeout',1e4);
set(ai,'SamplesPerTrigger',T2);
addchannel(ai,0:(Nchan-1));
ai.Channel.InputRange = [-10 10];
ai.Channel.SensorRange = [-10 10];

%ai.Channel.InputRange = [-5 5];
%ai.Channel.SensorRange = [-5 5];


%ai.Channel.SensorRange = [-10 10];

disp(['done opening ',filesave])

start(ai)


%data = zeros(T2,Nchan);
%data(:) = getdata(ai,T2);
clear('data');

data = getdata(ai,(T2 - 10*Fs));
%


%data = getdata(ai,ai.SamplesAvailable);


disp([' saving ',filesave])

save(filesave,'data','Fs');


disp(['done saving ',filesave])

%figure(gcf); clf;
%plot(data(1:2*Fs,3));

delete(ai);
clear a2

break


Ch_puff = 3;
Ch_stim = 2;
Ch_pupD = 1;
Ch_EEG = 4;

      

Tpre = 20;
Tpost = 20;

 t0 = [1:size(data,1)];

  
  Trigs = find(diff(data(:,Ch_stim)) > .5);
  Trigs(find(diff(Trigs)./Fs < 1)) = [];
  Trigs2  = Trigs./Fs;
  
  Npre = Tpre*Fs;
  Npost = Tpost*Fs;
  Npts = Npre + Npost;
  t = [[1:Npts]./Fs - Tpre]*1e3;
  
  Ntrials = length(Trigs);
  
  chans = 1:size(data,2); %[Ch_puff Ch_pupD Ch_EEG];
  Nchans = length(chans);
  Resp_Mat = zeros(Nchans, Ntrials, Npts);
  
  for count_trial = 1:Ntrials
    ind = find(t0 > (Trigs(count_trial) - Npre) & (t0 <= ...
						  Trigs(count_trial) + Npost));
    
    if length(ind) == Npts
      Resp_Mat(:,count_trial,:) = data(ind,chans)';
    end
    
  end
  
  
  
  
  
  

t2 = t./1e3; %./Fs;
ind_pre = find(t2 > -1 & t2<=0);
ind_post = find(t2>0.25 & t2 < .5);
ind_post2= find(t2>-.5 & t2 < .75);

Resp_pre = squeeze(mean(mean(Resp_Mat(Ch_pupD,:,ind_pre),1),3));
Resp_post = squeeze(mean(mean(Resp_Mat(Ch_pupD,:,ind_post),1),3));
Diff_Resp = (Resp_pre - Resp_post)./Resp_pre;

figure(gcf); clf;
subplot(2,2,1);

plot(Diff_Resp);

Ch_USE = Ch_pupD;
%Ch_USE = Ch_EEG;
%Ch_USE = Ch_stim;
%Ch_USE = Ch_puff;


ind_use2a = 2:2:Ntrials; %ind_use(find(Diff_Resp2>.25));

subplot(2,2,2)
tmp2a = squeeze(mean(Resp_Mat(Ch_USE,ind_use2a,:),2));
plot(t2(ind_post2),tmp2a(ind_post2),'b');

subplot(2,2,4)
tmp2a = squeeze(Resp_Mat(Ch_USE,ind_use2a,:));
plot(t2(ind_post2),tmp2a(:,ind_post2),'b');


break

