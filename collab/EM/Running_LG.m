%load and concatenate running data
data_base = 'G:\users\lindsey\dataLG';
DIR_base = fullfile(data_base,['LG' date '_' mouse], [date '_' mouse '\']);
run_data = [];
Trigs_all = [];
 Ntimestamps_tot = 0;
 Ntrigs_start_skip= 3;
     ch_TRIG = 6;
    ch_run = 3;

for iRun = 1:length(userun);
    DIR_data = dir([DIR_base, '*' [num2str(userun(iRun)) '.mat']]);
    load(fullfile(DIR_base,DIR_data.name));
    Ntimestamps0 = size(data,1);
    Trigs = find(diff(data(:,ch_TRIG))>1);
    Trigs(find(diff(Trigs)./Fs < .4)) = []; %in s
    Trigs(1: Ntrigs_start_skip) = [];   
    run_data = [run_data; data];
    Trigs_all = [Trigs_all; Trigs + Ntimestamps_tot];
    Ntimestamps_tot = Ntimestamps_tot + Ntimestamps0;
 end
clear('data');
%caclulate running speed during stim epochs
Nframes_pre = 10;
Nframes_post = 10;
Trigs_stim = Trigs_all( Nframes_pre+1 : (Nframes_pre + Nframes_post):end);
Nstim = length(Trigs_stim);


Trigs_run = find(run_data(:,ch_run)> .5);
Trigs_run(find(diff(Trigs_run)./Fs < .005));

WIN_RUN = [0 5];

SCALE_WHEEL = .048;

Wheel_mat = zeros(Nstim,1);
for count_stim = 1:Nstim
    ind = find(Trigs_run > (Trigs_stim(count_stim) - WIN_RUN(1)*Fs) & ...
        Trigs_run <= (Trigs_stim(count_stim) + WIN_RUN(2)*Fs));
    Wheel_mat(count_stim) = length(ind) ./ diff(WIN_RUN) *SCALE_WHEEL;
end

base_out = 'G:\users\lindsey\analysisLG\active mice';
fn_out = fullfile(base_out,mouse,date,[date '_' mouse '_run' num2str(userun) '_wheel.mat']);
save(fn_out, 'Wheel_mat');
clear('run_data');

%resort wheel data
seqfile = [date '_' mouse '_run' num2str(userun) '_Big_Seqposition.mat'];
load(fullfile(base_out,mouse,date,'analysis',seqfile));
load(fullfile(base_out,mouse,date,[date '_' mouse '_run' num2str(userun) '_wheel.mat']));

Wheel_sorted = zeros(size(Wheel_mat));
start = 1;
for iCond = 1:nCond;
    nRep = length(Big_Seqposition(iCond).ind);
    for iRep = 1:nRep;
        ind = Big_Seqposition(iCond).ind(iRep);
        Wheel_sorted(start) = Wheel_mat(ind);
        start = start+1;
    end
end

nblanks = length(Big_Seqposition(end).ind);
for iblank = 1:nblanks;
    ind = Big_Seqposition(end).ind(iblank);
    Wheel_sorted(start) = Wheel_mat(ind);
    start = start+1;
end

%resort stack data by running
%sort into running and nonrunning frames
stack_sorted = readtiff(fullfile(base_out,mouse,date,'analysis', [date '_' mouse '_run' num2str(userun) '_sorted.tif']));
load(fullfile(base_out,mouse,date,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']));
A = length(find(Wheel_mat>0));
B = length(find(Wheel_mat==0));

siz = size(stack_sorted);
stack_sorted_run = zeros(siz(1),siz(2),A*(nOFF+nON));
stack_sorted_norun = zeros(siz(1),siz(2),B*(nOFF+nON));

start_run = 0;
start_norun = 0;
start = 0;
stim_reps_run= zeros(1,nCond+1);
stim_reps_norun= zeros(1,nCond+1);
stim = 0;
for iCond = 1:nCond+1;
    count_run = 0;
    count_norun = 0;
    for itrial = 1+stim:stim_reps(iCond)+stim;
        if Wheel_mat(itrial)>0
            stack_sorted_run(:,:,1+start_run:start_run+nOFF+nON)= stack_sorted(:,:,1+start:start+nOFF+nON);
            start_run = start_run + nOFF+ nON;
            count_run = count_run+1;
        else
            stack_sorted_norun(:,:,1+start_norun:start_norun+nOFF+nON)= stack_sorted(:,:,1+start:start+nOFF+nON);
            start_norun = start_norun + nOFF+ nON;
            count_norun = count_norun+1;
        end
        start = start+nOFF+ nON;
    end
    stim = stim +stim_reps(iCond);
    stim_reps_run(1,iCond)= count_run;
    stim_reps_norun(1,iCond)= count_norun;
end

fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun) '_stim_reps_run.mat']);
save(fn_out, 'stim_reps_run');
fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun) '_stim_reps_norun.mat']);
save(fn_out, 'stim_reps_norun');
clear('stack_sorted');

stack_avg_run = zeros(siz(1), siz(2), nON+nOFF, nCond+1);
stack_avg_base_run = zeros(siz(1), siz(2), nCond+1);
stack_avg_resp_run = zeros(siz(1), siz(2), nCond+1);
stack_dF_run = zeros(siz(1), siz(2), nCond+1);
start = 0;
for iCond = 1:nCond+1;
    for iFrame = 1:(nOFF+nON)
    stack_avg_run(:,:,iFrame,iCond) = mean(stack_sorted_run(:,:,iFrame+start:nOFF+nON:start+((nOFF+nON)*stim_reps_run(iCond))),3);
    end
    stack_avg_base_run(:,:,iCond) = mean(stack_avg_run(:,:,pre_win(1):pre_win(2),iCond),3);
    stack_avg_resp_run(:,:,iCond) = mean(stack_avg_run(:,:,post_win(1):post_win(2),iCond),3);
    stack_dF_run(:,:,iCond) = (stack_avg_resp_run(:,:,iCond)-stack_avg_base_run(:,:,iCond))./stack_avg_base_run(:,:,iCond);
    start = (nOFF+nON)*stim_reps_run(iCond)+start;
end

stack_avg_norun = zeros(siz(1), siz(2), nON+nOFF, nCond+1);
stack_avg_base_norun = zeros(siz(1), siz(2), nCond+1);
stack_avg_resp_norun = zeros(siz(1), siz(2), nCond+1);
stack_dF_norun = zeros(siz(1), siz(2), nCond+1);
start = 0;
for iCond = 1:nCond+1;
    for iFrame = 1:(nOFF+nON)
    stack_avg_norun(:,:,iFrame,iCond) = mean(stack_sorted_norun(:,:,iFrame+start:nOFF+nON:start+((nOFF+nON)*stim_reps_norun(iCond))),3);
    end
    stack_avg_base_norun(:,:,iCond) = mean(stack_avg_norun(:,:,pre_win(1):pre_win(2),iCond),3);
    stack_avg_resp_norun(:,:,iCond) = mean(stack_avg_norun(:,:,post_win(1):post_win(2),iCond),3);
    stack_dF_norun(:,:,iCond) = (stack_avg_resp_norun(:,:,iCond)-stack_avg_base_norun(:,:,iCond))./stack_avg_base_norun(:,:,iCond);
    start = (nOFF+nON)*stim_reps_norun(iCond)+start;
end

fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun) '_stack_dF_run.tif']);
writetiff(stack_dF_run, fn_out);
fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun) '_stack_dF_norun.tif']);
writetiff(stack_dF_norun, fn_out);

clear('stack_sorted_run');
clear('stack_sorted_norun');
clear('stack_avg_run');
clear('stack_avg_norun');