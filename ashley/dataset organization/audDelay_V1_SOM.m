params.frameRate = 30;
params.stimOnTime = 100;
params.nBaselineMs = 1000;
params.nStimMs_visAlign = 3000;
params.nStimMs_audAlign = 9000;
params.nFramesVisDelay = 4;
params.motionCutoff = 0.2;

%% 931 180314
expt(1).SubNum = '931';
expt(1).mouse = '931';
expt(1).date = '180314';
expt(1).img_loc  = {'V1';'L2/3'};
expt(1).z = -284;
expt(1).img_strct  = {'cells'};
expt(1).indicator = {'virus';'SOM-GCaMP6s'};
expt(1).time_mat = ['1002';'1106'];
expt(1).runs = ['002';'004'];
expt(1).nrun = size(expt(1).runs,1);
expt(1).frame_rate = 30;
expt(1).rettuning = {'001';'0942'};
expt(1).stimOnMs = 100;
expt(1).motionThreshold = 0.05;
expt(1).areaBorders = 0;
expt(1).regImgStartFrame = 11800;
%% 931 180316
expt(2).SubNum = '931';
expt(2).mouse = '931';
expt(2).date = '180316';
expt(2).img_loc  = {'V1';'L2/3'};
expt(2).z = -284;
expt(2).img_strct  = {'cells'};
expt(2).indicator = {'virus';'SOM-GCaMP6s'};
expt(2).time_mat = ['1157';'1259'];
expt(2).runs = ['003';'004'];
expt(2).nrun = size(expt(2).runs,1);
expt(2).frame_rate = 30;
expt(2).rettuning = {'001';'1132'};
expt(2).stimOnMs = 100;
expt(2).motionThreshold = 0.05;
expt(2).areaBorders = 0;
expt(2).regImgStartFrame = '';
