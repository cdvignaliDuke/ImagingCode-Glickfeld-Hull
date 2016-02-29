%% AW13 150504
expt(1).SubNum = '613';
expt(1).mouse = 'AW13';
expt(1).date = '150504';
expt(1).img_loc  = {'V1';'L2/3'};
expt(1).img_strct  = {'cells'};
expt(1).time_mat = ['1425'; '1441';'1457'];
expt(1).runs = ['002'; '003'; '004'];
expt(1).nrun = size(expt(1).runs,1);
expt(1).frame_rate = 30;
expt(1).folder = 'two-photon imaging';
expt(1).catch = 0;
expt(1).dirtuning = '006';
expt(iexp).dirtuning_time = '1518';
expt(1).rettuning = {'005';'1514'};
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
expt(iexp).dirtuning_time = '1510';
expt(2).rettuning = {'005';'1506'};
expt(2).motionTestBimodal = 0;

%% AW13 150511 - FA rate high
expt(3).SubNum = '613';
expt(3).mouse = 'AW13';
expt(3).date = '150511';
expt(3).img_loc  = {'V1';'L2/3'};
expt(3).img_strct  = {'cells'};
expt(3).time_mat = ['1420'; '1437'; '1453'];
expt(3).runs = ['001'; '002'; '003'];
expt(3).nrun = size(expt(3).runs,1);
expt(3).frame_rate = 30;
expt(3).folder = 'two-photon imaging';
expt(3).catch = 0;
expt(3).dirtuning = '005';
expt(iexp).dirtuning_time = '1513';
expt(3).rettuning = {'004';'1509'};
expt(3).motionTestBimodal = 0;

%% AW14 150623
expt(4).SubNum = '614';
expt(4).mouse = 'AW14';
expt(4).date = '150623';
expt(4).img_loc  = {'V1';'L2/3'};
expt(4).img_strct  = {'cells'};
expt(4).time_mat = ['1142'; '1158'; '1217'];
expt(4).runs = ['002'; '003'; '004'];
expt(4).nrun = size(expt(4).runs,1);
expt(4).frame_rate = 30;
expt(4).folder = 'two-photon imaging';
expt(4).catch = 0;
expt(4).dirtuning = '006';
expt(iexp).dirtuning_time = '1239';
expt(4).rettuning = {'005';'1235'};
expt(4).motionTestBimodal = 1;

%% AW14 150626
expt(5).SubNum = '614';
expt(5).mouse = 'AW14';
expt(5).date = '150626';
expt(5).img_loc  = {'V1';'L2/3'};
expt(5).img_strct  = {'cells'};
expt(5).time_mat = ['1458'; '1514'; '1532'];
expt(5).runs = ['001'; '002'; '003'];
expt(5).nrun = size(expt(5).runs,1);
expt(5).frame_rate = 30;
expt(5).folder = 'two-photon imaging';
expt(5).catch = 0;
expt(5).dirtuning = '005';
expt(iexp).dirtuning_time = '1558';
expt(5).rettuning = {'004';'1553'};
expt(5).motionTestBimodal = 0;

%% AW14 150701 - might need to cut off trials for behavior
expt(6).SubNum = '614';
expt(6).mouse = 'AW14';
expt(6).date = '150701';
expt(6).img_loc  = {'V1';'L2/3'};
expt(6).img_strct  = {'cells'};
expt(6).time_mat = ['1556'; '1615'; '1633'];
expt(6).runs = ['002'; '003'; '004'];
expt(6).nrun = size(expt(6).runs,1);
expt(6).frame_rate = 30;
expt(6).folder = 'two-photon imaging';
expt(6).catch = 1;
expt(6).dirtuning = '006';
expt(iexp).dirtuning_time = '1659';
expt(6).rettuning = {'005';'1650'};
expt(6).motionTestBimodal = 0;

% %% AW16 151019 - miss rate high
% expt(7).SubNum = '616';
% expt(7).mouse = 'AW16';
% expt(7).date = '151019';
% expt(7).img_loc  = {'V1';'L2/3'};
% expt(7).img_strct  = {'cells'};
% expt(7).time_mat = ['1114'; '1131'; '1147'];
% expt(7).runs = ['003'; '004'; '005'];
% expt(7).nrun = size(expt(7).runs,1);
% expt(7).frame_rate = 30;
% expt(7).folder = 'two-photon imaging';
% expt(7).catch = 1;
% expt(7).dirtuning = '008';
% expt(7).rettuning = {'007';'1208'};
% expt(7).motionTestBimodal = 0;

% %% AW16 151028 - FA rate high
% expt(8).SubNum = '616';
% expt(8).mouse = 'AW16';
% expt(8).date = '151028';
% expt(8).img_loc  = {'V1';'L2/3'};
% expt(8).img_strct  = {'cells'};
% expt(8).time_mat = ['1041'; '1115'; '1130'];
% expt(8).runs = ['002'; '004'; '005'];
% expt(8).nrun = size(expt(8).runs,1);
% expt(8).frame_rate = 30;
% expt(8).folder = 'two-photon imaging';
% expt(8).catch = 1;
% expt(8).dirtuning = '007';
% expt(8).rettuning = {'006';'1146'};
% expt(8).motionTestBimodal = 0;
% 
% %% AW16 151030 - FA rate high
% expt(9).SubNum = '616';
% expt(9).mouse = 'AW16';
% expt(9).date = '151030';
% expt(9).img_loc  = {'V1';'L2/3'};
% expt(9).img_strct  = {'cells'};
% expt(9).time_mat = ['1050'; '1107'; '1123'];
% expt(9).runs = ['002'; '003'; '004'];
% expt(9).nrun = size(expt(9).runs,1);
% expt(9).frame_rate = 30;
% expt(9).folder = 'two-photon imaging';
% expt(9).catch = 1;
% expt(9).dirtuning = '006';
% expt(9).rettuning = {'005';'1140'};
% expt(9).motionTestBimodal = 1;

%% AW16 151105 - maybe suppressed, good behavior
expt(7).SubNum = '616';
expt(7).mouse = 'AW16';
expt(7).date = '151105';
expt(7).img_loc  = {'V1';'L2/3'};
expt(7).img_strct  = {'cells'};
expt(7).time_mat = ['1255'; '1311'; '1330'];
expt(7).runs = ['004'; '005'; '007'];
expt(7).nrun = size(expt(7).runs,1);
expt(7).frame_rate = 30;
expt(7).folder = 'two-photon imaging';
expt(7).catch = 1;
expt(7).dirtuning = '008';
expt(iexp).dirtuning_time = '1351';
expt(7).rettuning = {'008';'1347'};
expt(7).motionTestBimodal = 0;

% %% AW13 150430
% expt(8).SubNum = '613';
% expt(8).mouse = 'AW13';
% expt(8).date = '150430';
% expt(8).img_loc  = {'V1';'L2/3'};
% expt(8).img_strct  = {'cells'};
% expt(8).time_mat = ['1505'; '1521'; '1538'];
% expt(8).runs = [ '001'; '002'; '004'];
% expt(8).nrun = size(expt(8).runs,1);
% expt(8).frame_rate = 30;
% expt(8).folder = 'two-photon imaging';
% expt(8).catch = 0;
% expt(8).dirtuning = '006';
% expt(iexp).dirtuning_time = '1559';
% expt(8).rettuning = {'005';'1555'};
% expt(8).motionTestBimodal = 0;
% expt(8).z = 299;
% 
% %% AW13 150428
% expt(9).SubNum = '613';
% expt(9).mouse = 'AW13';
% expt(9).date = '150428';
% expt(9).img_loc  = {'V1';'L2/3'};
% expt(9).img_strct  = {'cells'};
% expt(9).time_mat = ['0955'; '1011'; '1029'];
% expt(9).runs = ['004'; '005'; '006'];
% expt(9).nrun = size(expt(9).runs,1);
% expt(9).frame_rate = 30;
% expt(9).folder = 'two-photon imaging';
% expt(9).catch = 0;
% expt(9).dirtuning = '008';
% expt(iexp).dirtuning_time = '1050';
% expt(9).rettuning = {'007';'1046'};
% expt(9).motionTestBimodal = 0;
% expt(9).z = 288;