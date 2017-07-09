%File info script for CRP sessions


dates = {'170330', '170416', '170426', '170427', '170417', '170425', ...
    '170426', '170427', '170420', '170428', '170427', '170428', '170429', '170429', '170501', '170501', '170505', '170506', '170501', ...
    '170509', '170508', '170510', '170513', '170515', '170518', '170519', '170519', '170520', '170522', '170523', '170524', ...
    '170529', '170530', '170531', '170606'};
mouseID = {'img86', 'img90', 'img90', 'img90', 'img91', 'img91', 'img91', 'img91', 'img92', ...
     'img92','img92', 'img90', 'img92', 'img90', 'img92', 'img90', 'img90', 'img90', 'img92', 'img92', 'img90', 'img93', 'img89', 'img90', ...
     'img93', 'img93', 'img89', 'img93', 'img89', 'img89', 'img94', 'img94', 'img94', 'img94', 'img94'};
 fileID = {};
behavID = {'989','990', '991', '992'};
runID = {'000', '001', '002'};
irun = 1;

%------------------------------------------------------
%day 1 of 500ms:  [2, 5, 9, 22, 23, 31]
days_1 = {'170416_img90',  '170417_img91', '170420_img92', '170510_img93', '170513_img89', '170524_img94'};

%post learning 500ms [3, 6, 10, 25, 27, 32]
days_post = {'170426_img90',  '170425_img91', '170428_img92', '170518_img93', '170519_img89', '170529_img94'};

%Unexpected Reward
days_UR = {'170506_img90', '170426_img91', '170429_img92', '170519_img93', '170522_img89', '170530_img94'};

%1000ms day 1
days_1000 = {'170508_img90', '170427_img91', '170501_img92', '170520_img93', '170523_img89', '170531_img94'};

%1000ms post learning
days_1000_post = {'170515_img90', [], [], [], [], '170606_img94'};

%----------------------------------------------------------
%img89
days89 = {'170513_img89', '170519_img89', '170522_img89', '170523_img89', []};

%img90
days90 = {'170416_img90', '170501_img90', '170506_img90', '170508_img90', '170515_img90'};  % '170426_img90'=10%omit   '170427_img90'=10%/10%   '170501_img90'=20%omit(500)  '170505_img90'=20%omit(500)quit early

%img91
days91 = {'170417_img91', '170425_img91', '170426_img91', '170427_img91', []}; %'170425_img91'=postLearning10%omit     '170426_img91'=10%/10%(500)

%img92
days92 = {'170420_img92', '170428_img92', '170429_img92', '170501_img92', []};   %'170427_img92'=10%omit(500) 

%img93
days93 = {'170510_img93', '170518_img93', '170519_img93', '170520_img93', []};

%img94
days94 = {'170524_img94', '170529_img94', '170530_img94', '170531_img94', '170606_img94'};

%%
             %89     90    91   92    93    94
comp_500ms = [23,27; 2:3; 5:6; 9,10; 22,25; 31,32];

%%
% day1/day_post      89              90             91             92             93             94
num_frames_to_use = [130000,130000; 130000,100000; 115000,120000; 120000,130000; 130000,130000; 130000,130000];



