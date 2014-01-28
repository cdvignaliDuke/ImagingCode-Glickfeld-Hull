%set params
RegInfo.All.path = '\\crash.dhe.duke.edu\lindsey\Data\2P_images\test\pollen_test';
RegInfo.All.fbase = {'pollen_test_run1'};
RegInfo.All.ch = '_green';
RegInfo.All.fnums = [1];
RegInfo.All.ftarget = 'pollen_test_target.tif';
RegInfo.All.is_rect = 0;
RegInfo.All.reverse_ch = 0;
RegInfo.All.nPlanes = 16;
Nvolbin = 1;
RegInfo.All.Nupsamp = 4;
RegInfo.All.Nfromedge = [6 6 6];
RegInfo.stack_fix_sz = {[512 512 400]};


matlabpool local6 % start the multiple virtual instances of matlab. "local6" is a 6 core local cluster profile I created using "Parallel" -> "Manage Cluster Profiles". 
trial_per_thread = 10;
RegInfo=doShift4D(RegInfo,trial_per_thread);
ShiftMat=RegInfo.All.ShiftMat;
figure;plot(ShiftMat(:,1))
figure;plot(ShiftMat(:,2))
figure;plot(ShiftMat(:,3))

%%
FrameChunk = 3000; % control memory usage
trial_per_thread = 10; % how many volumes to realign at a time in each core
RegInfo.All.limits=[-3 4 -3 8 0 1];
RegInfo.All.Nupsamp_XYZ=[1 1 2];
RegInfo=doReAlign4D(RegInfo,trial_per_thread,FrameChunk);