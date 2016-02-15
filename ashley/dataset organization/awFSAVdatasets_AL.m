%% AW26 160116
expt(1).SubNum = '626';
expt(1).mouse = 'AW26';
expt(1).date = '160116';
expt(1).img_loc  = {'AL';'L2/3'};
expt(1).img_strct  = {'cells'};
expt(1).time_mat = ['1536'; '1603'; '1618';'1637'];
expt(1).runs = ['001';'002'; '003'; '004'];
expt(1).nrun = size(expt(1).runs,1);
expt(1).frame_rate = 30;
expt(1).folder = 'two-photon imaging';
expt(1).catch = 1;
expt(1).dirtuning = '006';
expt(1).rettuning = {'005'; '1653'};
expt(1).motionTestBimodal = 0;

% %% AW26 160111 - questionable dataset, few cells, high FA rate
% expt(2).SubNum = '626';
% expt(2).mouse = 'AW26';
% expt(2).date = '160111';
% expt(2).img_loc  = {'AL';'L2/3'};
% expt(2).img_strct  = {'cells'};
% expt(2).time_mat = ['1422'; '1438'; '1455';'1511'];
% expt(2).runs = ['002'; '003'; '004';'005'];
% expt(2).nrun = size(expt(2).runs,1);
% expt(2).frame_rate = 30;
% expt(2).folder = 'two-photon imaging';
% expt(2).catch = 1;
% expt(2).dirtuning = '007';
% expt(2).rettuning = {'006'; '1528'};
% expt(2).motionTestBimodal = 0;

%% AW26 160203 - might be off-target location, responses somewhat suppressed overall, miss rate somewhat high
expt(2).SubNum = '626';
expt(2).mouse = 'AW26';
expt(2).date = '160203';
expt(2).img_loc  = {'AL';'L2/3'};
expt(2).img_strct  = {'cells'};
expt(2).time_mat = ['1053'; '1109'; '1125'];
expt(2).runs = ['001'; '002'; '003'];
expt(2).nrun = size(expt(2).runs,1);
expt(2).frame_rate = 30;
expt(2).folder = 'two-photon imaging';
expt(2).catch = 1;
expt(2).dirtuning = '005';
expt(2).rettuning = {'004'; '1142'};
expt(2).motionTestBimodal = 0;

%% AW26 160120 - very few cells, responses somewhat suppressed overall
expt(3).SubNum = '626';
expt(3).mouse = 'AW26';
expt(3).date = '160120';
expt(3).img_loc  = {'AL';'L2/3'};
expt(3).img_strct  = {'cells'};
expt(3).time_mat = ['1039'; '1057'; '1113';'1129'];
expt(3).runs = ['003'; '004'; '005';'006'];
expt(3).nrun = size(expt(3).runs,1);
expt(3).frame_rate = 30;
expt(3).folder = 'two-photon imaging';
expt(3).catch = 1;
expt(3).dirtuning = '010';
expt(3).rettuning = {'008'; '1147'};
expt(3).motionTestBimodal = 0;

%% AW26 160113 - very few cells, suppressesd by FSAV
expt(4).SubNum = '626';
expt(4).mouse = 'AW26';
expt(4).date = '160113';
expt(4).img_loc  = {'AL';'L2/3'};
expt(4).img_strct  = {'cells'};
expt(4).time_mat = ['1459'; '1515'; '1532';'1548'];
expt(4).runs = ['002'; '003'; '005';'006'];
expt(4).nrun = size(expt(4).runs,1);
expt(4).frame_rate = 30;
expt(4).folder = 'two-photon imaging';
expt(4).catch = 1;
expt(4).dirtuning = '008';
expt(4).rettuning = {'007'; '1604'};
expt(4).motionTestBimodal = 0;
