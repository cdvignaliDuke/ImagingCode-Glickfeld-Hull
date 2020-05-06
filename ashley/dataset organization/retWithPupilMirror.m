params.frameRate = 15.49;
params.nBaselineMs = 1000;
params.nPostStimMs = 2000;
params.nFramesVisDelay_VSR = 2;
params.motionCutoff = 0.2;
params.itiTimeS = 2;
params.eyeCalibMmPerPix = 1/26.6;
params.eyeCalibUmPerDeg = 25;
params.runSpeedThreshold_cps = 2;
params.eyePosThreshold_deg = 5;

%% 1103 190719
expt(1).mouse = '1103';
expt(1).date = '190719';
expt(1).img_loc  = {'PM';'L2/3'};
expt(1).z = -200;
expt(1).img_strct  = {'cells'};
expt(1).indicator = {'tg';''};
expt(1).retNoMirror = {'001','1613'};
expt(1).retWithMirror = {'003','1659'};
expt(1).mirrorPos = {'Az:0','El:0'};
expt(1).eyeradrange = [1 20];
expt(1).saveLoc = 'ashley';
expt(1).regImgStartFrame = 3721;

