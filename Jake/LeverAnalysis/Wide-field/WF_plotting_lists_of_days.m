%lists of WF experiments grouped according to various measures


%--------updated 2/26/17
%days = {'150518_img24', '150519_img24', '150518_img25', '150517_img25', '150716_img27', '150718_img27', '150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160516_img47'}; %'150718_img27', '150719_img27',
% days = {'150716_img28', '150717_img28', '151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46'}; %'150718_img27', '150719_img27',
% days = {'160129_img36','160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46'}; %'150718_img27', '150719_img27',
% %these are the days which have all their frames shrunk. Can use them to
% %analyze frame shift effect and look at licking. 
% days = {'160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46'}; %'150718_img27', '150719_img27',
% days = {'160921_img61', '160920_img61'};
days = {'161031_img68','161101_img68', '161030_img69', '161030_img70', '161101_img69', '161101_img70'};
%days = {'161108_img68', '161108_img69', '161108_img70'};  %cue-reward pairing delayed reward. 
days = {'161107_img68', '161107_img69', '161030_img68', '161030_img70', '161030_img69', '161031_img68', '161101_img69', '161101_img70'};  %first and last days

days = {'160209_img36', '151222_img32', '151019_img30', '160725_img53', '160905_img55'};  %no lever controls 161109_img61 161109_img59   Do not meet criterion 160208_img35

%days = {'160904_img55', '160905_img55', '160916_img61', '160918_img61', '160920_img61', '160921_img61', '161030_img62', '160904_img55'};

%updated 10/10/17
days = {'151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
days_chrono_order = {'151009_img30', '151011_img30','151021_img29', '151022_img29', '151211_img32', '151212_img32', '160129_img35', '160129_img36', '160131_img35', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
LS_ROIs =          {[2],               [2],           [1:3],           [1,2,3],        [1:2],          [1:3],        [1:3],           [1:3],           [1:3],         [1:3],           [1,2,3],         [1,2,3],        [1,2,3],        [1:2],           [3:4],        [1:4],          [1:5]}; %updated 10/10/17 LS
V_ROIs =          {[],                [],            [4,5],           [4,5],          [3:4],          [4:5],           [],              [],              [],           [],             [4],            [5],             [4],            [6],             [5],           [],             []}; %updated 10/10/17 V 
C1_ROIs =          {[1],               [1],           [],              [],               [],             [],           [],             [],              [],           [],             [5,6],            [4],            [5],            [4,5],           [],           [],               []}; %updated 10/10/17 C1
ROIcell = {        [2],               [2],           [1:3],            [1,3],           [1:4],           [1:5],       [1:2],           [1:2],           [1:2],         [1:2],           [3:4],           [2,3,5],         [1],            [1:2],         [3:4],        [1,2,3],         [3:5]}; % ROIs used for the main scatterplot in the manuscript

%Determined via by comparing the ratio of the peak corr df/f to that of the peak early df/f for only trials with licks. 
%Required to have a ratio of at least 1.15:1 corr:early. Calculated in
%scatter_lick_trials_only. Indeces refer to ROIs after removal of non-LS ROIs
valid_LS_ROIs =  {[],  [1],  [],  [],  [1,2],  [],  [1,2,3],  [],  [],  [2],  [],  [],  [], [], [1,2], [], [1,3,4]}; 

