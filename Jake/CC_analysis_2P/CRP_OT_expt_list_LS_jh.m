%CRP_OT_expt_list_LS_jh
%script containing the session info for all mice imaged on the overtraining
%paradigm in the lobule simplex

expt(1).name= 'Simp'; %day 1   
expt(1).mouse = strvcat('095', '096', '098', '1030', '1032');
expt(1).date = strvcat('190309', '190411', '190615', '190615', '190615');
expt(1).run = strvcat('000', '000', '000', '000', '000');
expt(1).ttl =   [1, 1, 1, 1, 1];
expt(1).noreg = [0, 0, 0, 0, 0];
expt(1).imgreglaseron = [0, 0, 0, 0, 0];

expt(2).name= 'Simp'; %OR PL  
expt(2).mouse = strvcat('095', '096', '098', '1032');
expt(2).date = strvcat('190313', '190415', '190621', '190621');
expt(2).run = strvcat('000', '000', '001', '000');
expt(2).ttl = [1, 1, 1, 1];
expt(2).noreg = [0, 0, 0, 0];
expt(2).imgreglaseron = [0, 0, 0, 0];

expt(3).name= 'Simp'; %UR   
expt(3).mouse = strvcat('095', '096');
expt(3).date = strvcat('190315', '190416');
expt(3).run = strvcat('000', '000');
expt(3).ttl = [1, 1];
expt(3).noreg = [0, 0];
expt(3).imgreglaseron = [0, 0];
%exp 5: trial ended right after targetOn (5 frames) on one trial

expt(4).name= 'Simp'; %Overtrained OR
expt(4).mouse =strvcat('095', '096', '098', '1030', '1032');
expt(4).date = strvcat('190321', '190423', '190628', '190628', '190628');
expt(4).run = strvcat('000', '000', '000', '000', '000');
expt(4).ttl = [1, 1, 1, 1, 1];
expt(4).noreg = [0, 0, 0, 0, 0];
expt(4).imgreglaseron = [0, 0, 0, 0, 0];

expt(5).name= 'Simp'; %Overtrained UR
expt(5).mouse =strvcat('095', '096');
expt(5).date = strvcat('190326', '190424');
expt(5).run = strvcat('000', '000');
expt(5).ttl = [1, 1];
expt(5).noreg = [0, 0];
expt(5).imgreglaseron = [0, 0];

%datasets w/o extractTC outputs

%datasets w/o spikeAlign

%datasets w/o lickAlign

%%datasets w/o lickResp

%datasets w/o lickResp_preCueOnly


