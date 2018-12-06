% %% AW16 151109 - ok (maybe best dataset), FA rate high
% expt(1).SubNum = '616';
% expt(1).mouse = 'AW16';
% expt(1).date = '151109';
% expt(1).img_loc  = {'V1';'L2/3'};
% expt(1).img_strct  = {'cells'};
% expt(1).time_mat = ['1102'; '1119'; '1137'];
% expt(1).runs = ['001'; '002'; '004'];
% expt(1).nrun = size(expt(1).runs,1);
% expt(1).frame_rate = 30;
% expt(1).folder = 'two-photon imaging';
% expt(1).catch = 1;
% expt(1).dirtuning = '007';
% expt(1).rettuning = {'005';'1154'};
% expt(1).motionTestBimodal = 1;
% expt(1).regImgStartFrame = 55416;
% 
% %% AW16 151111 - ok resp, few success trials, high FA rate
% expt(2).SubNum = '616';
% expt(2).mouse = 'AW16';
% expt(2).date = '151111';
% expt(2).img_loc  = {'V1';'L2/3'};
% expt(2).img_strct  = {'cells'};
% expt(2).time_mat = ['0851'; '0909'; '0925'];
% expt(2).runs = ['001'; '003'; '004'];
% expt(2).nrun = size(expt(2).runs,1);
% expt(2).frame_rate = 30;
% expt(2).folder = 'two-photon imaging';
% expt(2).catch = 1;
% expt(2).dirtuning = '006';
% expt(2).rettuning = {'005';'0948'};
% expt(2).motionTestBimodal = 0;
% expt(2).regImgStartFrame = 28710;

%% AW16 151124 - poor vis resp, but retinotopy should be good..., FA rate close to 50%
expt(1).SubNum = '616';
expt(1).mouse = 'AW16';
expt(1).date = '151124';
expt(1).img_loc  = {'V1';'L2/3'};
expt(1).img_strct  = {'cells'};
expt(1).time_mat = ['1047'; '1103'; '1119'];
expt(1).runs = ['001'; '002'; '003'];
expt(1).nrun = size(expt(1).runs,1);
expt(1).frame_rate = 30;
expt(1).folder = 'two-photon imaging';
expt(1).catch = 1;
expt(1).dirtuning = '005';
expt(1).dirtuning_time = '1145';
expt(1).rettuning = {'004';'1140'};
expt(1).motionTestBimodal = 0;
expt(1).trial_range = nan;
expt(1).regImgStartFrame = 41184;

% %% AW16 151116 - ok rsp, few trials, high FA rate
% expt(4).SubNum = '616';
% expt(4).mouse = 'AW16';
% expt(4).date = '151116';
% expt(4).img_loc  = {'V1';'L2/3'};
% expt(4).img_strct  = {'cells'};
% expt(4).time_mat = ['1044'; '1100'; '1116'];
% expt(4).runs = ['001'; '002'; '003'];
% expt(4).nrun = size(expt(4).runs,1);
% expt(4).frame_rate = 30;
% expt(4).folder = 'two-photon imaging';
% expt(4).catch = 1;
% expt(4).dirtuning = '005';
% expt(4).rettuning = {'004';'1133'};
% expt(4).motionTestBimodal = 0;
% expt(4).regImgStartFrame = 74437;
% 
% %% AW16 151118 - poor vis resp, high FA rate
% 
% expt(5).SubNum = '616';
% expt(5).mouse = 'AW16';
% expt(5).date = '151118';
% expt(5).img_loc  = {'V1';'L2/3'};
% expt(5).img_strct  = {'cells'};
% expt(5).time_mat = ['1048'; '1104'; '1122'];
% expt(5).runs = ['001'; '002'; '004'];
% expt(5).nrun = size(expt(5).runs,1);
% expt(5).frame_rate = 30;
% expt(5).folder = 'two-photon imaging';
% expt(5).catch = 1;
% expt(5).dirtuning = '006';
% expt(5).rettuning = {'005';'1138'};
% expt(5).motionTestBimodal = 0;
% expt(5).regImgStartFrame = 75825;

%% AW26 160112 - suppressed avg response, very few cells good bx
expt(2).SubNum = '626';
expt(2).mouse = 'AW26';
expt(2).date = '160112';
expt(2).img_loc  = {'V1';'L2/3'};
expt(2).img_strct  = {'cells'};
expt(2).time_mat = ['1453'; '1509'; '1525'; '1541'];
expt(2).runs = ['001'; '002'; '003'; '004'];
expt(2).nrun = size(expt(2).runs,1);
expt(2).frame_rate = 30;
expt(2).folder = 'two-photon imaging';
expt(2).catch = 1;
expt(2).dirtuning = '006';
expt(2).dirtuning_time = '1602';
expt(2).rettuning = {'005'; '1558'};
expt(2).motionTestBimodal = 0;
expt(2).trial_range = 50:390;
expt(2).regImgStartFrame = 67956;
%% AW26 160118 - suppressed avg response, very few cells, good behav
expt(3).SubNum = '626';
expt(3).mouse = 'AW26';
expt(3).date = '160118';
expt(3).img_loc  = {'V1';'L2/3'};
expt(3).img_strct  = {'cells'};
expt(3).time_mat = ['1617'; '1650'; '1707'];
expt(3).runs = ['003'; '004'; '005'];
expt(3).nrun = size(expt(3).runs,1);
expt(3).frame_rate = 30;
expt(3).folder = 'two-photon imaging';
expt(3).catch = 1;
expt(3).dirtuning = '007';
expt(3).dirtuning_time = '1729';
expt(3).rettuning = {'006'; '1724'};
expt(3).motionTestBimodal = 0;
expt(3).trial_range = 1:250;
expt(3).regImgStartFrame = 58624;
% 
% %% AW26 160125 - suppressed avg response, very few cells, good behavior
% expt(8).SubNum = '626';
% expt(8).mouse = 'AW26';
% expt(8).date = '160125';
% expt(8).img_loc  = {'V1';'L2/3'};
% expt(8).img_strct  = {'cells'};
% expt(8).time_mat = ['1416'; '1432'; '1449'; '1505'];
% expt(8).runs = ['001'; '002'; '003'; '004'];
% expt(8).nrun = size(expt(8).runs,1);
% expt(8).frame_rate = 30;
% expt(8).folder = 'two-photon imaging';
% expt(8).catch = 1;
% expt(8).dirtuning = '006';
% expt(8).rettuning = {'005';'1522'};
% expt(8).motionTestBimodal = 0;
% expt(8).regImgStartFrame = 1587;
% 
% %% AW26 160127 - suppressed avg response, very few cells, high miss rate
% expt(9).SubNum = '626';
% expt(9).mouse = 'AW26';
% expt(9).date = '160127';
% expt(9).img_loc  = {'V1';'L2/3'};
% expt(9).img_strct  = {'cells'};
% expt(9).time_mat = ['1116'; '1135'; '1151'; '1210'];
% expt(9).runs = ['001'; '002'; '003'; '004'];
% expt(9).nrun = size(expt(9).runs,1);
% expt(9).frame_rate = 30;
% expt(9).folder = 'two-photon imaging';
% expt(9).catch = 1;
% expt(9).dirtuning = '006';
% expt(9).rettuning = {'005';'1228'};
% expt(9).motionTestBimodal = 0;
% expt(9).regImgStartFrame = '';

%% AW26 160202 **BEST AW26 V1 dataset - suppressed avg response, very few cells, miss rate somewhat high (could cut beginning trials)
expt(4).SubNum = '626';
expt(4).mouse = 'AW26';
expt(4).date = '160202';
expt(4).img_loc  = {'V1';'L2/3'};
expt(4).img_strct  = {'cells'};
expt(4).time_mat = ['1335'; '1350'; '1407'];
expt(4).runs = ['003'; '004'; '006'];
expt(4).nrun = size(expt(4).runs,1);
expt(4).frame_rate = 30;
expt(4).folder = 'two-photon imaging';
expt(4).catch = 1;
expt(4).dirtuning = '008';
expt(4).dirtuning_time = '1435';
expt(4).rettuning = {'007';'1424'};
expt(4).motionTestBimodal = 0;
expt(4).trial_range = 74:290;
expt(4).regImgStartFrame = 69612;

% %% AW26 160121 - suppressed avg response, no cells, good behavior
% expt(5).SubNum = '626';
% expt(5).mouse = 'AW26';
% expt(5).date = '160121';
% expt(5).img_loc  = {'V1';'L2/3'};
% expt(5).img_strct  = {'cells'};
% expt(5).time_mat = ['1045'; '1100'; '1117'];
% expt(5).runs = ['001'; '002'; '003'];
% expt(5).nrun = size(expt(5).runs,1);
% expt(5).frame_rate = 30;
% expt(5).folder = 'two-photon imaging';
% expt(5).catch = 1;
% expt(5).dirtuning = '006';
% expt(5).dirtuning_time = '1215';
% expt(5).rettuning = {'005';'1202'};
% expt(5).motionTestBimodal = 0;
% expt(5).trial_range = 1:250;
% expt(5).regImgStartFrame = 58624;

%% AW25 160302 - high FA rate
expt(5).SubNum = '625';
expt(5).mouse = 'AW25';
expt(5).date = '160302';
expt(5).img_loc  = {'V1';'L2/3'};
expt(5).img_strct  = {'cells'};
expt(5).time_mat = ['1113';'1130'];%'1040'; '1057';
expt(5).runs = ['003';'004'];%'001'; '002';
expt(5).nrun = size(expt(5).runs,1);
expt(5).frame_rate = 30;
expt(5).folder = 'two-photon imaging';
expt(5).catch = 1;
expt(5).dirtuning = '006';
expt(5).dirtuning_time = '1150';
expt(5).rettuning = {'005';'1146'};
expt(5).motionTestBimodal = 0;
expt(5).z = 281;
expt(5).trial_range = nan;
expt(5).regImgStartFrame = 29136;
%% AW25 160309
expt(6).SubNum = '625';
expt(6).mouse = 'AW25';
expt(6).date = '160309';
expt(6).img_loc  = {'V1';'L2/3'};
expt(6).img_strct  = {'cells'};
expt(6).time_mat = ['1346';'1402';'1421';'1437'];%
expt(6).runs = ['006';'007';'008';'009'];% 
expt(6).nrun = size(expt(6).runs,1);
expt(6).frame_rate = 30;
expt(6).folder = 'two-photon imaging';
expt(6).catch = 1;
expt(6).dirtuning = '011';
expt(6).dirtuning_time = '1457';
expt(6).rettuning = {'010';'1453'};
expt(6).motionTestBimodal = 0;
expt(6).z = 244;
expt(6).trial_range = nan;
expt(6).regImgStartFrame = 94528;

% %% AW25 150311 - FA rate somewhat high
% expt(14).SubNum = '625';
% expt(14).mouse = 'AW25';
% expt(14).date = '160311';
% expt(14).img_loc  = {'V1';'L2/3'};
% expt(14).img_strct  = {'cells'};
% expt(14).time_mat = ['1026'; '1042';'1058'];%'1005'; 
% expt(14).runs = ['002'; '003';'004'];%'001'; 
% expt(14).nrun = size(expt(14).runs,1);
% expt(14).frame_rate = 30;
% expt(14).folder = 'two-photon imaging';
% expt(14).catch = 1;
% expt(14).dirtuning = '006';
% expt(14).dirtuning_time = '1121';
% expt(14).rettuning = {'005';'1116'};
% expt(14).motionTestBimodal = 0;
% expt(14).z = 237;
% expt(14).trial_range = 26:230;
% expt(14).regImgStartFrame = 14141;

%% AW25 150314 **run
expt(7).SubNum = '625';
expt(7).mouse = 'AW25';
expt(7).date = '160314';
expt(7).img_loc  = {'V1';'L2/3'};
expt(7).img_strct  = {'cells'};
expt(7).time_mat = ['1021'; '1037'];%;'1053'];
expt(7).runs = [ '002'; '003'];%;'004'];
expt(7).nrun = size(expt(7).runs,1);
expt(7).frame_rate = 30;
expt(7).folder = 'two-photon imaging';
expt(7).catch = 1;
expt(7).dirtuning = '007';
expt(7).dirtuning_time = '1119';
expt(7).rettuning = {'006';'1115'};
expt(7).motionTestBimodal = 0;
expt(7).z = 245;
expt(7).trial_range = nan;
expt(7).regImgStartFrame = 50742;

%% AW25 160316
expt(8).SubNum = '625';
expt(8).mouse = 'AW25';
expt(8).date = '160316';
expt(8).img_loc  = {'V1';'L2/3'};
expt(8).img_strct  = {'cells'};
expt(8).time_mat = ['1224'; '1240';'1256'];%;'1312'];
expt(8).runs = ['001';'002'; '003'];%; '004'];
expt(8).nrun = size(expt(8).runs,1);
expt(8).frame_rate = 30;
expt(8).folder = 'two-photon imaging';
expt(8).catch = 0;
expt(8).catchRew = 0;
expt(8).dirtuning = '006';
expt(8).dirtuning_time = '1333';
expt(8).rettuning = {'005';'1328'};
expt(8).z = 267;
expt(8).regImgStartFrame = 3253;
expt(8).trial_range = 1:250;
expt(8).regImgStartFrame = 26280;

% %% AW16 151019 - miss rate high
% expt(17).SubNum = '616';
% expt(17).mouse = 'AW16';
% expt(17).date = '151019';
% expt(17).img_loc  = {'V1';'L2/3'};
% expt(17).img_strct  = {'cells'};
% expt(17).time_mat = ['1114'; '1131'; '1147'];
% expt(17).runs = ['003'; '004'; '005'];
% expt(17).nrun = size(expt(17).runs,1);
% expt(17).frame_rate = 30;
% expt(17).folder = 'two-photon imaging';
% expt(17).catch = 1;
% expt(17).dirtuning = '008';
% expt(17).rettuning = {'007';'1208'};
% expt(17).motionTestBimodal = 0;
% expt(17).regImgStartFrame = '';
% 
% %% AW16 151028 - FA rate high
% expt(18).SubNum = '616';
% expt(18).mouse = 'AW16';
% expt(18).date = '151028';
% expt(18).img_loc  = {'V1';'L2/3'};
% expt(18).img_strct  = {'cells'};
% expt(18).time_mat = ['1041'; '1115'; '1130'];
% expt(18).runs = ['002'; '004'; '005'];
% expt(18).nrun = size(expt(18).runs,1);
% expt(18).frame_rate = 30;
% expt(18).folder = 'two-photon imaging';
% expt(18).catch = 1;
% expt(18).dirtuning = '007';
% expt(18).rettuning = {'006';'1146'};
% expt(18).motionTestBimodal = 0;
% 
% %% AW16 151030 - FA rate high
% expt(19).SubNum = '616';
% expt(19).mouse = 'AW16';
% expt(19).date = '151030';
% expt(19).img_loc  = {'V1';'L2/3'};
% expt(19).img_strct  = {'cells'};
% expt(19).time_mat = ['1050'; '1107'; '1123'];
% expt(19).runs = ['002'; '003'; '004'];
% expt(19).nrun = size(expt(19).runs,1);
% expt(19).frame_rate = 30;
% expt(19).folder = 'two-photon imaging';
% expt(19).catch = 1;
% expt(19).dirtuning = '006';
% expt(19).rettuning = {'005';'1140'};
% expt(19).motionTestBimodal = 1;

%% AW16 151105 - maybe suppressed, good behavior
expt(9).SubNum = '616';
expt(9).mouse = 'AW16';
expt(9).date = '151105';
expt(9).img_loc  = {'V1';'L2/3'};
expt(9).img_strct  = {'cells'};
expt(9).time_mat = ['1255'; '1311'; '1330'];
expt(9).runs = ['004'; '005'; '007'];
expt(9).nrun = size(expt(9).runs,1);
expt(9).frame_rate = 30;
expt(9).folder = 'two-photon imaging';
expt(9).catch = 1;
expt(9).catchRew = 1;
expt(9).dirtuning = '008';
expt(9).dirtuning_time = '1351';
expt(9).rettuning = {'008';'1347'};
expt(9).motionTestBimodal = 0;
expt(9).trial_range = [30:140,171:270];
expt(9).regImgStartFrame = 54988;