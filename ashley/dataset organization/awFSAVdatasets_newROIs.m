%% AW14 150623
expt(1).SubNum = '614';
expt(1).mouse = 'AW14';
expt(1).date = '150623';
expt(1).img_loc  = {'V1';'L2/3'};
expt(1).img_strct  = {'cells'};
expt(1).time_mat = ['1142'; '1158'; '1217'];
expt(1).runs = ['002'; '003'; '004'];
expt(1).nrun = size(expt(1).runs,1);
expt(1).frame_rate = 30;
expt(1).folder = 'two-photon imaging';
expt(1).catch = 0;
expt(1).dirtuning = '006';
expt(1).dirtuning_time = '1239';
expt(1).rettuning = {'005';'1235'};
expt(1).motionTestBimodal = 1;

%% AW13 150508
expt(2).SubNum = '613';
expt(2).mouse = 'AW13';
expt(2).date = '150508';
expt(2).img_loc  = {'V1';'L2/3'};
expt(2).img_strct  = {'cells'};
expt(2).time_mat = ['1413'; '1430'; '1446'];
expt(2).runs = ['001'; '002'; '003'];
expt(2).nrun = size(expt(2).runs,1);
expt(2).frame_rate = 30;
expt(2).folder = 'two-photon imaging';
expt(2).catch = 1;
expt(2).dirtuning = '006';
expt(2).dirtuning_time = '1510';
expt(2).rettuning = {'005';'1506'};
expt(2).motionTestBimodal = 0;

%% AW25 160316 - 100ms
expt(3).SubNum = '625';
expt(3).mouse = 'AW25';
expt(3).date = '160316';
expt(3).img_loc  = {'V1';'L2/3'};
expt(3).img_strct  = {'cells'};
expt(3).time_mat = ['1224'; '1240'; '1256';'1312'];
expt(3).runs = ['001'; '002'; '003';'004'];
expt(3).nrun = size(expt(3).runs,1);
expt(3).frame_rate = 30;
expt(3).folder = 'two-photon imaging';
expt(3).catch = 1;
expt(3).dirtuning = '006';
expt(3).dirtuning_time = '1333';
expt(3).rettuning = {'005';'1328'};
expt(3).motionTestBimodal = 0;
expt(3).z = 267;