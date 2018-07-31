LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

%%
date = '170923';
ImgFolder = strvcat('001','002','003');
time = strvcat('1247','1319','1350');
mouse = 'i698';
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% did not use- eye partly closed often
date = '170923';
ImgFolder = strvcat('001');
time = '1440';
mouse = 'i699';
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% 
date = '171013';
ImgFolder = strvcat('001');
time = '1011';
mouse = 'i699';
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%%
date = '170127';
ImgFolder = strvcat('002','003','004');
time = strvcat('1152','1211','1233');
mouse = 'i671';
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);
%% did not use
date = '170131';
ImgFolder = strvcat('003','004','005','006');
time = strvcat('1445','1502','1517','1553');
mouse = 'i671';
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);
