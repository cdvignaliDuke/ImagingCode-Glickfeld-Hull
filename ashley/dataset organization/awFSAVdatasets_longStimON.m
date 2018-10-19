%% AW16 151109 - ok (maybe best dataset), FA rate high
expt(1).SubNum = '616';
expt(1).mouse = 'AW16';
expt(1).date = '151109';
expt(1).img_loc  = {'V1';'L2/3'};
expt(1).img_strct  = {'cells'};
expt(1).time_mat = ['1102'; '1119'; '1137'];
expt(1).runs = ['001'; '002'; '004'];
expt(1).nrun = size(expt(1).runs,1);
expt(1).frame_rate = 30;
expt(1).folder = 'two-photon imaging';
expt(1).catch = 1;
expt(1).dirtuning = '007';
expt(1).rettuning = {'005';'1154'};
expt(1).motionTestBimodal = 1;

%% AW16 151111 - ok resp, few success trials, high FA rate
expt(2).SubNum = '616';
expt(2).mouse = 'AW16';
expt(2).date = '151111';
expt(2).img_loc  = {'V1';'L2/3'};
expt(2).img_strct  = {'cells'};
expt(2).time_mat = ['0851'; '0909'; '0925'];
expt(2).runs = ['001'; '003'; '004'];
expt(2).nrun = size(expt(2).runs,1);
expt(2).frame_rate = 30;
expt(2).folder = 'two-photon imaging';
expt(2).catch = 1;
expt(2).dirtuning = '006';
expt(2).rettuning = {'005';'0948'};
expt(2).motionTestBimodal = 0;

%% AW16 151124 - poor vis resp, but retinotopy should be good..., FA rate close to 50%
expt(3).SubNum = '616';
expt(3).mouse = 'AW16';
expt(3).date = '151124';
expt(3).img_loc  = {'V1';'L2/3'};
expt(3).img_strct  = {'cells'};
expt(3).time_mat = ['1047'; '1103'; '1119'];
expt(3).runs = ['001'; '002'; '003'];
expt(3).nrun = size(expt(3).runs,1);
expt(3).frame_rate = 30;
expt(3).folder = 'two-photon imaging';
expt(3).catch = 1;
expt(3).dirtuning = '005';
expt(3).dirtuning_time = '1145';
expt(3).rettuning = {'004';'1140'};
expt(3).motionTestBimodal = 0;
expt(3).trial_range = nan;

%% AW16 151116 - ok rsp, few trials, high FA rate
expt(4).SubNum = '616';
expt(4).mouse = 'AW16';
expt(4).date = '151116';
expt(4).img_loc  = {'V1';'L2/3'};
expt(4).img_strct  = {'cells'};
expt(4).time_mat = ['1044'; '1100'; '1116'];
expt(4).runs = ['001'; '002'; '003'];
expt(4).nrun = size(expt(4).runs,1);
expt(4).frame_rate = 30;
expt(4).folder = 'two-photon imaging';
expt(4).catch = 1;
expt(4).dirtuning = '005';
expt(4).rettuning = {'004';'1133'};
expt(4).motionTestBimodal = 0;

%% AW16 151118 - poor vis resp, high FA rate

expt(5).SubNum = '616';
expt(5).mouse = 'AW16';
expt(5).date = '151118';
expt(5).img_loc  = {'V1';'L2/3'};
expt(5).img_strct  = {'cells'};
expt(5).time_mat = ['1048'; '1104'; '1122'];
expt(5).runs = ['001'; '002'; '004'];
expt(5).nrun = size(expt(5).runs,1);
expt(5).frame_rate = 30;
expt(5).folder = 'two-photon imaging';
expt(5).catch = 1;
expt(5).dirtuning = '006';
expt(5).rettuning = {'005';'1138'};
expt(5).motionTestBimodal = 0;

%% AW26 160112 - suppressed avg response, very few cells good bx
expt(6).SubNum = '626';
expt(6).mouse = 'AW26';
expt(6).date = '160112';
expt(6).img_loc  = {'V1';'L2/3'};
expt(6).img_strct  = {'cells'};
expt(6).time_mat = ['1453'; '1509'; '1525'; '1541'];
expt(6).runs = ['001'; '002'; '003'; '004'];
expt(6).nrun = size(expt(6).runs,1);
expt(6).frame_rate = 30;
expt(6).folder = 'two-photon imaging';
expt(6).catch = 1;
expt(6).dirtuning = '006';
expt(6).dirtuning_time = '1602';
expt(6).rettuning = {'005'; '1558'};
expt(6).motionTestBimodal = 0;
expt(6).trial_range = nan;
%% AW26 160118 - suppressed avg response, very few cells, good behav
expt(7).SubNum = '626';
expt(7).mouse = 'AW26';
expt(7).date = '160118';
expt(7).img_loc  = {'V1';'L2/3'};
expt(7).img_strct  = {'cells'};
expt(7).time_mat = ['1617'; '1650'; '1707'];
expt(7).runs = ['003'; '004'; '005'];
expt(7).nrun = size(expt(7).runs,1);
expt(7).frame_rate = 30;
expt(7).folder = 'two-photon imaging';
expt(7).catch = 1;
expt(7).dirtuning = '007';
expt(7).dirtuning_time = '1729';
expt(7).rettuning = {'006'; '1724'};
expt(7).motionTestBimodal = 0;
expt(7).trial_range = nan;

%% AW26 160125 - suppressed avg response, very few cells, good behavior
expt(8).SubNum = '626';
expt(8).mouse = 'AW26';
expt(8).date = '160125';
expt(8).img_loc  = {'V1';'L2/3'};
expt(8).img_strct  = {'cells'};
expt(8).time_mat = ['1416'; '1432'; '1449'; '1505'];
expt(8).runs = ['001'; '002'; '003'; '004'];
expt(8).nrun = size(expt(8).runs,1);
expt(8).frame_rate = 30;
expt(8).folder = 'two-photon imaging';
expt(8).catch = 1;
expt(8).dirtuning = '006';
expt(8).rettuning = {'005';'1522'};
expt(8).motionTestBimodal = 0;

%% AW26 160127 - suppressed avg response, very few cells, high miss rate
expt(9).SubNum = '626';
expt(9).mouse = 'AW26';
expt(9).date = '160127';
expt(9).img_loc  = {'V1';'L2/3'};
expt(9).img_strct  = {'cells'};
expt(9).time_mat = ['1116'; '1135'; '1151'; '1210'];
expt(9).runs = ['001'; '002'; '003'; '004'];
expt(9).nrun = size(expt(9).runs,1);
expt(9).frame_rate = 30;
expt(9).folder = 'two-photon imaging';
expt(9).catch = 1;
expt(9).dirtuning = '006';
expt(9).rettuning = {'005';'1228'};
expt(9).motionTestBimodal = 0;

%% AW26 160202 **BEST AW26 V1 dataset - suppressed avg response, very few cells, miss rate somewhat high (could cut beginning trials)
expt(10).SubNum = '626';
expt(10).mouse = 'AW26';
expt(10).date = '160202';
expt(10).img_loc  = {'V1';'L2/3'};
expt(10).img_strct  = {'cells'};
expt(10).time_mat = ['1335'; '1350'; '1407'];
expt(10).runs = ['003'; '004'; '006'];
expt(10).nrun = size(expt(10).runs,1);
expt(10).frame_rate = 30;
expt(10).folder = 'two-photon imaging';
expt(10).catch = 1;
expt(10).dirtuning = '008';
expt(10).dirtuning_time = '1435';
expt(10).rettuning = {'007';'1424'};
expt(10).motionTestBimodal = 0;
expt(10).trial_range = nan;

%% AW26 160121 - suppressed avg response, no cells, good behavior
expt(11).SubNum = '626';
expt(11).mouse = 'AW26';
expt(11).date = '160121';
expt(11).img_loc  = {'V1';'L2/3'};
expt(11).img_strct  = {'cells'};
expt(11).time_mat = ['1045'; '1100'; '1117'; '1133'];
expt(11).runs = ['001'; '002'; '003'; '004'];
expt(11).nrun = size(expt(11).runs,1);
expt(11).frame_rate = 30;
expt(11).folder = 'two-photon imaging';
expt(11).catch = 1;
expt(11).dirtuning = '006';
expt(11).rettuning = {'005';'1202'};
expt(11).motionTestBimodal = 0;

%% AW25 160302 - high FA rate
expt(12).SubNum = '625';
expt(12).mouse = 'AW25';
expt(12).date = '160302';
expt(12).img_loc  = {'V1';'L2/3'};
expt(12).img_strct  = {'cells'};
expt(12).time_mat = ['1113';'1130'];%'1040'; '1057';
expt(12).runs = ['003';'004'];%'001'; '002';
expt(12).nrun = size(expt(12).runs,1);
expt(12).frame_rate = 30;
expt(12).folder = 'two-photon imaging';
expt(12).catch = 1;
expt(12).dirtuning = '006';
expt(12).dirtuning_time = '1150';
expt(12).rettuning = {'005';'1146'};
expt(12).motionTestBimodal = 0;
expt(12).z = 281;
expt(12).trial_range = nan;
%% AW25 160309
expt(13).SubNum = '625';
expt(13).mouse = 'AW25';
expt(13).date = '160309';
expt(13).img_loc  = {'V1';'L2/3'};
expt(13).img_strct  = {'cells'};
expt(13).time_mat = ['1346';'1402';'1421';'1437'];%
expt(13).runs = ['006';'007';'008';'009'];% 
expt(13).nrun = size(expt(13).runs,1);
expt(13).frame_rate = 30;
expt(13).folder = 'two-photon imaging';
expt(13).catch = 1;
expt(13).dirtuning = '011';
expt(13).dirtuning_time = '1457';
expt(13).rettuning = {'010';'1453'};
expt(13).motionTestBimodal = 0;
expt(13).z = 244;
expt(13).trial_range = nan;

%% AW25 150311 - FA rate somewhat high
expt(14).SubNum = '625';
expt(14).mouse = 'AW25';
expt(14).date = '160311';
expt(14).img_loc  = {'V1';'L2/3'};
expt(14).img_strct  = {'cells'};
expt(14).time_mat = ['1026'; '1042';'1058'];%'1005'; 
expt(14).runs = ['002'; '003';'004'];%'001'; 
expt(14).nrun = size(expt(14).runs,1);
expt(14).frame_rate = 30;
expt(14).folder = 'two-photon imaging';
expt(14).catch = 1;
expt(14).dirtuning = '006';
expt(14).dirtuning_time = '1121';
expt(14).rettuning = {'005';'1116'};
expt(14).motionTestBimodal = 0;
expt(14).z = 237;
expt(14).trial_range = nan;

%% AW25 150314 **run
expt(15).SubNum = '625';
expt(15).mouse = 'AW25';
expt(15).date = '160314';
expt(15).img_loc  = {'V1';'L2/3'};
expt(15).img_strct  = {'cells'};
expt(15).time_mat = ['1005'; '1021'; '1037';'1053'];
expt(15).runs = ['001'; '002'; '003';'004'];
expt(15).nrun = size(expt(15).runs,1);
expt(15).frame_rate = 30;
expt(15).folder = 'two-photon imaging';
expt(15).catch = 1;
expt(15).dirtuning = '007';
expt(15).dirtuning_time = '1119';
expt(15).rettuning = {'006';'1115'};
expt(15).motionTestBimodal = 0;
expt(15).z = 245;
