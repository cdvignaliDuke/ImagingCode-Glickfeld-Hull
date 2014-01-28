plane_number = 1; %1 => 1st plane (or volume) imaged that day

PWD = 'G:\users\lindsey\analysisLG\active mice\AC42\110928\';
PWD_EEG = 'G:\users\lindsey\dataLG\2010\110210\LG25\';

file_EEG_mat = strvcat(...
    'eegacq_mouse_training2011-2-9_16-11-40_rec_110209_CFJ_run1.mat', ...
    'eegacq_mouse_training2011-2-9_16-49-43_rec_110209_CFJ_run2.mat', ...
    'eegacq_mouse_training2011-2-9_17-19-20_rec_110209_CFJ_run3.mat', ...
    'eegacq_mouse_training2011-2-9_17-55-15_rec_110209_CFJ_run4.mat', ...
    'eegacq_mouse_training2011-2-9_18-24-37_rec_110209_CFJ_run5.mat');

file_mat = strvcat(...
    '110928_AC42_run1.tif', ...
    '110928_AC42_run2.tif', ...
    '110928_AC42_run3.tif', ...
    '110928_AC42_run4.tif', ...
    '110928_AC42_run5.tif', ...
    '110928_AC42_run6.tif', ...
    '110928_AC42_run7.tif', ...
    '110928_AC42_run8.tif');

%file_mat = ['100502_u3912_run3023__ch2_fix_reg.tif'];
file_ch2 = 'AVG_110210_LG25_run1.tif';
file_mask = 'AVG_110210_LG25_run1.tif';

DIR_rand0 = strvcat('G:\users\lindsey\dataLG\vis_stim_protocols\110928\AC42_run2\TF5_SF5_ori1_Con1_Nfrpstim240_Ret1_GabSz20_28Sep2011\',...
    'G:\users\lindsey\dataLG\vis_stim_protocols\110928\AC42_run6\TF3_SF3_ori1_Con1_Nfrpstim240_Ret1_GabSz20_28Sep2011\');
                
Nprotocols = size(DIR_rand0,1);

randvec = [0 0
           1 1
           1 2
           1 3
           1 4
           2 1
           2 2
           2 3]; %directory #, randnum (within that directory)

RandStim_info = [];

%%%Repeat entry of these parameters for each stimulus protocol presented (i.e. Nprotocols): 
count_prot = 1;
RandStim_info.prot(count_prot).userun0 = [2:5];
RandStim_info.prot(count_prot).userunDR = [1];

RandStim_info.prot(count_prot).TFvec = [repmat([1 2 4 8 15],1,5) 0];
RandStim_info.prot(count_prot).SFvec = [.32*ones(1,5) .16*ones(1,5) .08*ones(1,5) .04*ones(1,5) .02*ones(1,5) 0];
RandStim_info.prot(count_prot).pos = [ones(1,26)];
RandStim_info.prot(count_prot).StimXpos = [0*ones(1,25) 0]; %x pos 
RandStim_info.prot(count_prot).StimYpos = [-15*ones(1,25) 0]; %y pos 
RandStim_info.prot(count_prot).Nstimtypes0000 = 26; %2 or 6 or 8
RandStim_info.prot(count_prot).oris = [90*ones(1,25) 0]*pi/180;
RandStim_info.prot(count_prot).contrastvec = [.8*ones(1,25) 0];
RandStim_info.prot(count_prot).Fs1 = 2; %acq rate of imaging Hz %downsampled if applicable
RandStim_info.prot(count_prot).Fs1_orig = 2; %acq rate of imaging Hz; ignore unless you have eeg
RandStim_info.prot(count_prot).Noff1 = 10;
RandStim_info.prot(count_prot).Non1 = 10;
RandStim_info.prot(count_prot).Toff1 = RandStim_info.prot(count_prot).Noff1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ton1 = RandStim_info.prot(count_prot).Non1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ntot1 = RandStim_info.prot(count_prot).Fs1 * ...
    (RandStim_info.prot(count_prot).Ton1 + RandStim_info.prot(count_prot).Toff1);


%%% in case of second protocol on the same day, fill in the following: 
if Nprotocols > 1
    count_prot = 2;
RandStim_info.prot(count_prot).userun0 = [6:8];
RandStim_info.prot(count_prot).userunDR = [1];

RandStim_info.prot(count_prot).TFvec = [repmat([1 4 15],1,3) 0];
RandStim_info.prot(count_prot).SFvec = [.32*ones(1,3) .08*ones(1,3) .02*ones(1,3) 0];
RandStim_info.prot(count_prot).pos = [ones(1,10)];
RandStim_info.prot(count_prot).StimXpos = [-30*ones(1,9) 0]; %x pos 
RandStim_info.prot(count_prot).StimYpos = [-15*ones(1,9) 0]; %y pos 
RandStim_info.prot(count_prot).Nstimtypes0000 = 10; %2 or 6 or 8
RandStim_info.prot(count_prot).oris = [90*ones(1,9) 0]*pi/180;
RandStim_info.prot(count_prot).contrastvec = [.8*ones(1,9) 0];
RandStim_info.prot(count_prot).Fs1 = 2; %acq rate of imaging Hz %downsampled if applicable
RandStim_info.prot(count_prot).Fs1_orig = 2; %acq rate of imaging Hz; ignore unless you have eeg
RandStim_info.prot(count_prot).Noff1 = 10;
RandStim_info.prot(count_prot).Non1 = 10;
RandStim_info.prot(count_prot).Toff1 = RandStim_info.prot(count_prot).Noff1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ton1 = RandStim_info.prot(count_prot).Non1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ntot1 = RandStim_info.prot(count_prot).Fs1 * ...
    (RandStim_info.prot(count_prot).Ton1 + RandStim_info.prot(count_prot).Toff1);end

if Nprotocols > 2
    count_prot = 3;
RandStim_info.prot(count_prot).userun0 = [7:9];
RandStim_info.prot(count_prot).userunDR = [1];

RandStim_info.prot(count_prot).TFvec = [repmat([1 1 2 2 4 4 8 8 15 15],1,5) 0];
RandStim_info.prot(count_prot).SFvec = [.32*ones(1,10) .16*ones(1,10) .08*ones(1,10) .04*ones(1,10) .02*ones(1,10) 0];
RandStim_info.prot(count_prot).pos = [ones(1,51)];
RandStim_info.prot(count_prot).StimXpos = [25*ones(1,50) 0]; %x pos 
RandStim_info.prot(count_prot).StimYpos = [0*ones(1,50) 0]; %y pos 
RandStim_info.prot(count_prot).Nstimtypes0000 = 51; %2 or 6 or 8
RandStim_info.prot(count_prot).oris = [repmat([90 270],1,25) 0]*pi/180;
RandStim_info.prot(count_prot).contrastvec = [.8*ones(1,50) 0];
RandStim_info.prot(count_prot).Fs1 = 2.667; %acq rate of imaging Hz %downsampled if applicable
RandStim_info.prot(count_prot).Fs1_orig = 2.667; %acq rate of imaging Hz; ignore unless you have eeg
RandStim_info.prot(count_prot).Noff1 = 12;
RandStim_info.prot(count_prot).Non1 = 12;
RandStim_info.prot(count_prot).Toff1 = RandStim_info.prot(count_prot).Noff1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ton1 = RandStim_info.prot(count_prot).Non1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ntot1 = RandStim_info.prot(count_prot).Fs1 * ...
    (RandStim_info.prot(count_prot).Ton1 + RandStim_info.prot(count_prot).Toff1);
end

if Nprotocols > 3
    count_prot = 4;
RandStim_info.prot(count_prot).userun0 = [6];
RandStim_info.prot(count_prot).userunDR = [1];

RandStim_info.prot(count_prot).TFvec = [1*ones(1,3) 0];
RandStim_info.prot(count_prot).SFvec = [0.16*ones(1,3) 0];
RandStim_info.prot(count_prot).pos = [ones(1,4)];
RandStim_info.prot(count_prot).StimXpos = [-20 0 20 0]; %x pos 
RandStim_info.prot(count_prot).StimYpos = [-10*ones(1,3) 0]; %y pos 
RandStim_info.prot(count_prot).Nstimtypes0000 =4; %2 or 6 or 8
RandStim_info.prot(count_prot).oris = [90*ones(1,3) 0]*pi/180;
RandStim_info.prot(count_prot).contrastvec = [.8*ones(1,3) 0];
RandStim_info.prot(count_prot).Fs1 = 2; %acq rate of imaging Hz %downsampled if applicable
RandStim_info.prot(count_prot).Fs1_orig = 2; %acq rate of imaging Hz; ignore unless you have eeg
RandStim_info.prot(count_prot).Noff1 = 10;
RandStim_info.prot(count_prot).Non1 = 10;
RandStim_info.prot(count_prot).Toff1 = RandStim_info.prot(count_prot).Noff1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ton1 = RandStim_info.prot(count_prot).Non1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ntot1 = RandStim_info.prot(count_prot).Fs1 * ...
    (RandStim_info.prot(count_prot).Ton1 + RandStim_info.prot(count_prot).Toff1);
end

if Nprotocols > 4
    count_prot = 5;
RandStim_info.prot(count_prot).userun0 = [7:9];
RandStim_info.prot(count_prot).userunDR = [1];
RandStim_info.prot(count_prot).TFvec = [repmat([1 1 2 2 4 4 8 8 15 15],1,5) 0];
RandStim_info.prot(count_prot).SFvec = [.32*ones(1,10) .16*ones(1,10) .08*ones(1,10) .04*ones(1,10) .02*ones(1,10) 0];
RandStim_info.prot(count_prot).pos = [ones(1,51)];
RandStim_info.prot(count_prot).StimXpos = [-10*ones(1,50) 0]; %x pos 
RandStim_info.prot(count_prot).StimYpos = [-10*ones(1,50) 0]; %y pos 
RandStim_info.prot(count_prot).Nstimtypes0000 = 51; %2 or 6 or 8
RandStim_info.prot(count_prot).oris = [repmat([90 270],1,25) 0]*pi/180;
RandStim_info.prot(count_prot).contrastvec = [.8*ones(1,50) 0];
RandStim_info.prot(count_prot).Fs1 = 2.667; %acq rate of imaging Hz %downsampled if applicable
RandStim_info.prot(count_prot).Fs1_orig = 2.667; %acq rate of imaging Hz; ignore unless you have eeg
RandStim_info.prot(count_prot).Noff1 = 12;
RandStim_info.prot(count_prot).Non1 = 12;
RandStim_info.prot(count_prot).Toff1 = RandStim_info.prot(count_prot).Noff1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ton1 = RandStim_info.prot(count_prot).Non1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ntot1 = RandStim_info.prot(count_prot).Fs1 * ...
    (RandStim_info.prot(count_prot).Ton1 + RandStim_info.prot(count_prot).Toff1);
end

if Nprotocols > 5
    count_prot = 6;
RandStim_info.prot(count_prot).userun0 = [10:12];
RandStim_info.prot(count_prot).userunDR = [1];
RandStim_info.prot(count_prot).TFvec = [repmat([1 1 2 2 4 4 8 8 15 15],1,5) 0];
RandStim_info.prot(count_prot).SFvec = [.32*ones(1,10) .16*ones(1,10) .08*ones(1,10) .04*ones(1,10) .02*ones(1,10) 0];
RandStim_info.prot(count_prot).pos = [ones(1,51)];
RandStim_info.prot(count_prot).StimXpos = [-20*ones(1,50) 0]; %x pos 
RandStim_info.prot(count_prot).StimYpos = [-15*ones(1,50) 0]; %y pos 
RandStim_info.prot(count_prot).Nstimtypes0000 = 51; %2 or 6 or 8
RandStim_info.prot(count_prot).oris = [repmat([90 270],1,25) 0]*pi/180;
RandStim_info.prot(count_prot).contrastvec = [.8*ones(1,50) 0];
RandStim_info.prot(count_prot).Fs1 = 2.667; %acq rate of imaging Hz %downsampled if applicable
RandStim_info.prot(count_prot).Fs1_orig = 2.667; %acq rate of imaging Hz; ignore unless you have eeg
RandStim_info.prot(count_prot).Noff1 = 12;
RandStim_info.prot(count_prot).Non1 = 12;
RandStim_info.prot(count_prot).Toff1 = RandStim_info.prot(count_prot).Noff1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ton1 = RandStim_info.prot(count_prot).Non1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ntot1 = RandStim_info.prot(count_prot).Fs1 * ...
    (RandStim_info.prot(count_prot).Ton1 + RandStim_info.prot(count_prot).Toff1);
end
%%%%% next, enter all the info above so it's easy to access later if needed
RandStim_info.general.PWD = PWD;
RandStim_info.general.PWD_EEG = PWD_EEG;
RandStim_info.general.file_EEG_mat = file_EEG_mat;
RandStim_info.general.file_mat = file_mat;
RandStim_info.general.file_ch2 = file_ch2;
RandStim_info.general.file_mask = file_mask;
RandStim_info.general.DIR_rand0 = DIR_rand0;
RandStim_info.general.Nprotocols = size(DIR_rand0,1);
RandStim_info.general.randvec = randvec;
RandStim_info.general.plane_number = plane_number; 
save([PWD,'RandStim_info_plane',num2str(plane_number),'.mat'],'RandStim_info');

