plane_number = 1; %1 => 1st plane (or volume) imaged that day

PWD = 'G:\users\lindsey\analysisLG\active mice\Y2\110508\';
PWD_EEG = 'G:\users\lindsey\dataLG\2010\110210\LG25\';

file_EEG_mat = strvcat(...
    'eegacq_mouse_training2011-2-9_16-11-40_rec_110209_CFJ_run1.mat', ...
    'eegacq_mouse_training2011-2-9_16-49-43_rec_110209_CFJ_run2.mat', ...
    'eegacq_mouse_training2011-2-9_17-19-20_rec_110209_CFJ_run3.mat', ...
    'eegacq_mouse_training2011-2-9_17-55-15_rec_110209_CFJ_run4.mat', ...
    'eegacq_mouse_training2011-2-9_18-24-37_rec_110209_CFJ_run5.mat');

file_mat = strvcat(...
    '110508_Y2_run1_max_reg.tif', ...
    '110508_Y2_run2_max_reg.tif', ...
    '110508_Y2_run3_max_reg.tif', ...
    '110508_Y2_run4_max_reg.tif', ...
    '110508_Y2_run5_max_reg.tif', ...
    '110508_Y2_run6_max_reg.tif');

%file_mat = ['100502_u12_run3023__ch2_fix_reg.tif'];
file_ch2 = 'AVG_110210_LG25_run1.tif';
file_mask = 'AVG_110210_LG25_run1.tif';

DIR_rand0 = strvcat('G:\users\lindsey\dataLG\vis_stim_protocols\110508\Y2_run1\TF6_SF6_ori1_Con1_Nfrpstim240_Ret1_GabSz13_08May2011');
                
Nprotocols = size(DIR_rand0,1);

randvec = [1 1
           1 2 
           1 3
           0 0
           2 1
           2 2]; %directory #, randnum (within that directory)
          

RandStim_info = [];

%%%Repeat entry of these parameters for each stimulus protocol presented (i.e. Nprotocols): 
count_prot = 1;
RandStim_info.prot(count_prot).userun0 = [1];
RandStim_info.prot(count_prot).userunDR = [1];

RandStim_info.prot(count_prot).TFvec = [2*ones(1,24) 8*ones(1,24) 0];
RandStim_info.prot(count_prot).SFvec = [repmat([.04*ones(1,8) .08*ones(1,8) .16*ones(1,8)],1,2) 0];
RandStim_info.prot(count_prot).pos = [ones(1,49)];
RandStim_info.prot(count_prot).StimXpos = [-5*ones(1,48) 0]; %x pos 
RandStim_info.prot(count_prot).StimYpos = [-15*ones(1,48) 0]; %y pos 
RandStim_info.prot(count_prot).Nstimtypes0000 = 49; %2 or 6 or 8
RandStim_info.prot(count_prot).oris = [repmat([0 45 90 135 180 225 270 315],1,6) 0]*pi/180;
RandStim_info.prot(count_prot).contrastvec = [.8*ones(1,48) 0];
RandStim_info.prot(count_prot).Fs1 = 4; %acq rate of imaging Hz %downsampled if applicable
RandStim_info.prot(count_prot).Fs1_orig = 4; %acq rate of imaging Hz; ignore unless you have eeg
RandStim_info.prot(count_prot).Noff1 = 18;
RandStim_info.prot(count_prot).Non1 = 18;
RandStim_info.prot(count_prot).Toff1 = RandStim_info.prot(count_prot).Noff1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ton1 = RandStim_info.prot(count_prot).Non1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ntot1 = RandStim_info.prot(count_prot).Fs1 * ...
    (RandStim_info.prot(count_prot).Ton1 + RandStim_info.prot(count_prot).Toff1);


%%% in case of second protocol on the same day, fill in the following: 
if Nprotocols > 1
    count_prot = 2;
RandStim_info.prot(count_prot).userun0 = [5];
RandStim_info.prot(count_prot).userunDR = [1];

RandStim_info.prot(count_prot).TFvec = [2 8 2*ones(1,8) 8*ones(1,8) 0];
RandStim_info.prot(count_prot).SFvec = [0.02 0.02 0.32*ones(1,16) 0];
RandStim_info.prot(count_prot).pos = [ones(1,19)];
RandStim_info.prot(count_prot).StimXpos = [-5*ones(1,18) 0]; %x pos 
RandStim_info.prot(count_prot).StimYpos = [-15*ones(1,18) 0]; %y pos 
RandStim_info.prot(count_prot).Nstimtypes0000 =19; %2 or 6 or 8
RandStim_info.prot(count_prot).oris = [0 0 0 45 90 135 180 225 270 315 0 45 90 135 180 225 270 315 0]*pi/180;
RandStim_info.prot(count_prot).contrastvec = [.8*ones(1,18) 0];
RandStim_info.prot(count_prot).Fs1 = 4; %acq rate of imaging Hz %downsampled if applicable
RandStim_info.prot(count_prot).Fs1_orig = 4; %acq rate of imaging Hz; ignore unless you have eeg
RandStim_info.prot(count_prot).Noff1 = 18;
RandStim_info.prot(count_prot).Non1 = 18;
RandStim_info.prot(count_prot).Toff1 = RandStim_info.prot(count_prot).Noff1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ton1 = RandStim_info.prot(count_prot).Non1/RandStim_info.prot(count_prot).Fs1; %s
RandStim_info.prot(count_prot).Ntot1 = RandStim_info.prot(count_prot).Fs1 * ...
    (RandStim_info.prot(count_prot).Ton1 + RandStim_info.prot(count_prot).Toff1);
end

if Nprotocols > 2
    count_prot = 3;
RandStim_info.prot(count_prot).userun0 = [6];
RandStim_info.prot(count_prot).userunDR = [1];

RandStim_info.prot(count_prot).TFvec = [8*ones(1,5) 0];
RandStim_info.prot(count_prot).SFvec = [0.02 0.04 0.08 0.16 0.32 0];
RandStim_info.prot(count_prot).pos = [ones(1,6)];
RandStim_info.prot(count_prot).StimXpos = [-10*ones(1,5) 0]; %x pos 
RandStim_info.prot(count_prot).StimYpos = [0*ones(1,5) 0]; %y pos 
RandStim_info.prot(count_prot).Nstimtypes0000 = 6; %2 or 6 or 8
RandStim_info.prot(count_prot).oris = [0*ones(1,5) 0]*pi/180;
RandStim_info.prot(count_prot).contrastvec = [.8*ones(1,5) 0];
RandStim_info.prot(count_prot).Fs1 = 4; %acq rate of imaging Hz %downsampled if applicable
RandStim_info.prot(count_prot).Fs1_orig = 4; %acq rate of imaging Hz; ignore unless you have eeg
RandStim_info.prot(count_prot).Noff1 = 18;
RandStim_info.prot(count_prot).Non1 = 18;
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

