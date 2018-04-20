
load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
%% V1 cell bodies, FS single interval, during behavior
date = '150626';
ImgFolder = strvcat('001','002','003');
time = strvcat('1458','1514', '1532');
mouse = 'i614';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% V1 cell bodies, FS single interval, during behavior
date = '161020';
ImgFolder = strvcat('001','002','003','004');
time = strvcat('1244','1300', '1315', '1332');
mouse = 'i669';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% V1 cell bodies, FS single interval, during behavior
date = '170116';
ImgFolder = strvcat('003','004','005');
time = strvcat('1114', '1131','1147');
mouse = 'i671';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% V1 cell bodies, FS random interval, passive
date = '170112';
ImgFolder = strvcat('002','003');
time = strvcat('1710','1734');
mouse = 'i674';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% V1 cell bodies, FS random interval, passive
date = '170116';
ImgFolder = strvcat('003','004');
time = strvcat('1347','1405');
mouse = 'i674';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% V1 axons in lateral area, FS random interval, passive
date = '170124';
ImgFolder = strvcat('002','003','004');
time = strvcat('1551','1608','1626');
mouse = 'i671';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% V1 cell bodies, FS random interval, during behavior
date = '170127';
ImgFolder = strvcat('002','003','004');
time = strvcat('1152','1211','1233');
mouse = 'i671';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% V1 cell bodies, FS random interval, during behavior
date = '170131';
ImgFolder = strvcat('003','004','005','006');
time = strvcat('1445','1502','1517','1533');
mouse = 'i671';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% V1 axons in medial area, FS random interval, passive
date = '170201';
ImgFolder = strvcat('004','005','006');
time = strvcat('1457','1515','1532');
mouse = 'i671';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% V1 cell bodies, FS random interval, passive
date = '170207';
ImgFolder = strvcat('004','005');
time = strvcat('1519','1536');
mouse = 'i684';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% V1 cell bodies, FS random interval, passive
date = '170210';
ImgFolder = strvcat('002','003');
time = strvcat('1626','1643');
mouse = 'i696';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% LM cell bodies, FS random interval, passive
date = '170220';
ImgFolder = strvcat('004');
time = strvcat('1717');
mouse = 'i684';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% PM cell bodies, FS random interval, passive
date = '170301';
ImgFolder = strvcat('002');
time = strvcat('1704');
mouse = 'i684';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);

%% PM cell bodies, FS random interval, passive
date = '170314';
ImgFolder = strvcat('002');
time = strvcat('1748');
mouse = 'i696';
doFromRef = 0;
ref = strvcat('005');
nrun = size(ImgFolder,1);
frame_rate = 30;
run_str = catRunName(ImgFolder, nrun);
