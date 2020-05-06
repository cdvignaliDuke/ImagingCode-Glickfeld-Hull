% 16 Direction across days- SLC mice
frame_rate = 15.5;
fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';

%% i1312
expt(1).mouse = 'i1312';
expt(1).ref_date = '200118';
expt(1).ref_run = strvcat('002');
expt(1).reg_dates = strvcat('200120','200201');
expt(1).reg_runs = strvcat('003','002');

expt(1).ndate = size(expt(1).reg_dates,1);

%% i1313
expt(2).mouse = 'i1313';
expt(2).ref_date = '200118';
expt(2).ref_run = strvcat('002');
expt(2).reg_dates = strvcat('200120','200201');
expt(2).reg_runs = strvcat('002','002');

expt(2).ndate = size(expt(2).reg_dates,1);
