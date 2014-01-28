mouse = 'AC45';
date = '110823';
userun = 1:4;

nCond = 50;

nON = 12;
nOFF = 12;
pre_win = [7 12];
post_win = [13 24];

base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse, date);

%fastrig running
DIR_base = fullfile('\\zoloto\bigstorlab\fastrig\running\Running_data',[date '_' mouse '\']);
ch_wheel = 1;

Fs_fastrig = 30.82; %Hz
dTrig0 = 1/Fs_fastrig;

Nframes_pre = 144;
Nframes_post = 144;
Wheel_mat = [];
eval(['PARAMS_' date '_' mouse]);
data_base = 'G:\users\lindsey\analysisLG\active mice';

for iRun  = 1:length(userun);
%figure out #frames
file_USE = deblank(file_mat(userun(iRun),:));
Nframes1 = sizetiff(fullfile(data_base,mouse,date,file_USE))*12;               

%load running data
DIR_data = dir([DIR_base, '*' [num2str(userun(iRun)) '.mat']]);
load(fullfile(DIR_base,DIR_data.name));
Ntimestamps0 = size(data,1);

%ID trigs
Trigs = ceil(Fs*[dTrig0:dTrig0:(dTrig0*Nframes1)]);

%Find eye mvmt and wheel movement in the same period
thresh_wheelmvmt = 3;
Trigs_wheel = find(diff(data(:,ch_wheel)) > thresh_wheelmvmt);
Trigs_wheel(find(diff(Trigs_wheel)./Fs < .01)) = [];

SCALE_WHEEL = .048; %10 ticks per revolution, 6" diam -> .48m circumference, so .048 m/tick

Nstim = Nframes1./(Nframes_pre+Nframes_post);
Wheel_mat_temp = zeros(Nstim,1);
for count = 1:Nstim-1
   ind = find(Trigs_wheel>(Trigs(1+Nframes_pre+((count-1)*(Nframes_pre+Nframes_post)))) & Trigs_wheel<=(Trigs(count*(Nframes_pre+Nframes_post))));
   Wheel_mat_temp(count) = length(ind)./(Nframes_post/Fs_fastrig) * SCALE_WHEEL;
end

Wheel_mat = [Wheel_mat; Wheel_mat_temp];
end

fn_out = fullfile(outDir,'analysis', [date '_' mouse '_run' num2str(userun) '_wheel.mat']);
save(fn_out, 'Wheel_mat');
clear('run_data');

%resort wheel data
seqfile = [date '_' mouse '_run' num2str(userun) '_Big_Seqposition.mat'];
load(fullfile(outDir,'analysis',seqfile));
load(fullfile(outDir,'analysis', [date '_' mouse '_run' num2str(userun) '_wheel.mat']));

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
stack_sorted = readtiff(fullfile(base,mouse,date,'analysis', [date '_' mouse '_run' num2str(userun) '_sorted.tif']));
load(fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']));
thresh_run = 0;
A = length(find(Wheel_mat>thresh_run));
B = length(find(Wheel_mat==thresh_run));

siz = size(stack_sorted);
stack_sorted_run = zeros(siz(1),siz(2),A*(nOFF+nON));
stack_sorted_norun = zeros(siz(1),siz(2),B*(nOFF+nON));

start_run = 0;
start_norun = 0;
start = 0;
stim_reps_run= zeros(1,nCond+1);
stim_reps_norun= zeros(1,nCond+1);
stim = 0;
run_ind = [];
norun_ind = [];
for iCond = 1:nCond+1;
    count_run = 0;
    count_norun = 0;
    for itrial = 1+stim:stim_reps(iCond)+stim;
        if Wheel_sorted(itrial)>thresh_run
            stack_sorted_run(:,:,1+start_run:start_run+nOFF+nON)= stack_sorted(:,:,1+start:start+nOFF+nON);
            start_run = start_run + nOFF+ nON;
            count_run = count_run+1;
            run_ind = [run_ind itrial];
        else
            stack_sorted_norun(:,:,1+start_norun:start_norun+nOFF+nON)= stack_sorted(:,:,1+start:start+nOFF+nON);
            start_norun = start_norun + nOFF+ nON;
            count_norun = count_norun+1;
            norun_ind = [norun_ind itrial];
        end
        start = start+nOFF+ nON;
    end
    stim = stim +stim_reps(iCond);
    stim_reps_run(1,iCond)= count_run;
    stim_reps_norun(1,iCond)= count_norun;
end

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_stim_reps_run.mat']);
save(fn_out, 'stim_reps_run');
fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_stim_reps_norun.mat']);
save(fn_out, 'stim_reps_norun');
fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_running_ind.mat']);
save(fn_out, 'run_ind', 'norun_ind');

clear('stack_sorted');

siz = size(stack_sorted_run);
stack_dF_run = zeros(siz(1), siz(2), nCond+1);
start = 0;
for iCond = 1:nCond+1;
    nRep = stim_reps_run(iCond);
    rep_dF = zeros(siz(1), siz(2),nRep);
    for iRep = 1:nRep
        rep_base = mean(stack_sorted_run(:,:,start+pre_win(1):start+pre_win(2)),3);
        rep_resp = mean(stack_sorted_run(:,:,start+post_win(1):start+post_win(2)),3);
        rep_dF(:,:,iRep) = (rep_resp-rep_base)./rep_base;
        start = nOFF+nON+start;
    end
    stack_dF_run(:,:,iCond) = mean(rep_dF,3);
end


siz = size(stack_sorted_norun);
stack_dF_norun = zeros(siz(1), siz(2), nCond+1);
start = 0;
for iCond = 1:nCond+1;
    nRep = stim_reps_norun(iCond);
    rep_dF = zeros(siz(1), siz(2),nRep);
    for iRep = 1:nRep
        rep_base = mean(stack_sorted_norun(:,:,start+pre_win(1):start+pre_win(2)),3);
        rep_resp = mean(stack_sorted_norun(:,:,start+post_win(1):start+post_win(2)),3);
        rep_dF(:,:,iRep) = (rep_resp-rep_base)./rep_base;
        start = nOFF+nON+start;
    end
    stack_dF_norun(:,:,iCond) = mean(rep_dF,3);
end

fn_out = fullfile(outDir,'analysis', [date '_' mouse '_run' num2str(userun) '_stack_dF_run.tif']);
writetiff(stack_dF_run, fn_out);
fn_out = fullfile(outDir, 'analysis',  [date '_' mouse '_run' num2str(userun) '_stack_dF_norun.tif']);
writetiff(stack_dF_norun, fn_out);

clear('stack_sorted_run');
clear('stack_sorted_norun');
clear('stack_avg_run');
clear('stack_avg_norun');

%%

Here is a chunk where I find all running onsets (and offsets or ends of runningbouts) where there is no running in the prv 2s, and at least some running in the next 2 s:


WIN_ONSET = [-3 8];
t_onset = WIN_ONSET(1):1/Fs1:WIN_ONSET(2);
t_onset(end) = [];

WIN_ONSET_PRE_SOMERUN = 2; %s
WIN_ONSET_POST_SOMERUN = 2; %something during each of 3 sec pre

WIN_OFFSET = WIN_ONSET
WIN_OFFSET_PRE_SOMERUN = 2; %s
WIN_OFFSET_POST_SOMERUN = 2; %s

ind_run = find(tmp1_wheel);
ind_run_onset = [];%ind_run;
ind_run_offset = [];%ind_run;

for count = 1:length(ind_run)
   ind_run0 = ind_run(count);
   ind_post = ind_run0 + [0:(WIN_ONSET_POST_SOMERUN*Fs1-1)];
   ind_pre = ind_run0 + [-1*(WIN_ONSET_PRE_SOMERUN*Fs1):-1];

   if min(ind_pre)>0 & min(ind_post)>0 & max(ind_pre)<=length(tmp1_wheel) & max(ind_post)<=length(tmp1_wheel)
       tmp1_wheel_post = sum(reshape(tmp1_wheel(ind_post),Fs1,WIN_ONSET_POST_SOMERUN),1)>0;
       tmp1_wheel_pre = sum(reshape(tmp1_wheel(ind_pre),Fs1,WIN_ONSET_PRE_SOMERUN),1)>0;
       if (sum(tmp1_wheel_post) == WIN_ONSET_POST_SOMERUN & sum(tmp1_wheel_pre) == 0)
           ind_run_onset = [ind_run_onset;  ind_run0];


       end
   end

   %repeat for offset
      ind_run0 = ind_run(count);
   ind_post = ind_run0 + [1:(WIN_OFFSET_POST_SOMERUN*Fs1)];
   ind_pre = ind_run0 +  1 + [-1*(WIN_OFFSET_PRE_SOMERUN*Fs1):-1];

   if min(ind_pre)>0 & min(ind_post)>0 & max(ind_pre)<=length(tmp1_wheel) & max(ind_post)<=length(tmp1_wheel)
       tmp1_wheel_post = sum(reshape(tmp1_wheel(ind_post),Fs1,WIN_OFFSET_POST_SOMERUN),1)>0;
       tmp1_wheel_pre = sum(reshape(tmp1_wheel(ind_pre),Fs1,WIN_OFFSET_PRE_SOMERUN),1)>0;
       if (sum(tmp1_wheel_pre) == WIN_OFFSET_PRE_SOMERUN & sum(tmp1_wheel_post) == 0)
           ind_run_offset = [ind_run_offset;  ind_run0];
       end
   end

end
%%%%%%%%%%%%%%%%%%%%%%


%can combine with mask so you look only at moments with or without vis stim..

Mask_isstim = zeros(length(tmp1_wheel),1);
Mask_isstim = reshape(Mask_isstim,Ntot1,length(tmp1_wheel)./Ntot1);
Mask_isstim(Noff1+1:end,:) = 1;
%find blank trials:
%Mask_isstim(:,find(Seq2(:,4)==0)) == 0;

Mask_isstim = reshape(Mask_isstim,length(tmp1_wheel),1);


tmp1_onset0 = (tmp1_wheel>0 & [1; tmp1_wheel(1:end-1)]==0 & [ones(2,1); tmp1_wheel(1:end-2)]==0 & [ones(3,1); tmp1_wheel(1:end-3)]==0);
tmp1_offset0 = ([1; tmp1_wheel(1:end-1)]>0 & [tmp1_wheel]==0 & [tmp1_wheel(2:end); 1]==0 & [tmp1_wheel(3:end); ones(2,1)]==0);
Mask_isstim_meanlum_pre3_post3 = (Mask_isstim==0 & [Mask_isstim(2:end); 1]==0 & [1; Mask_isstim(1:end-1)]==0);
Mask_isstim_meanlum_pre3_post3 = (Mask_isstim==0 & [Mask_isstim(2:end); 1]==0 & [Mask_isstim(3:end); ones(2,1)]==0 & [1; Mask_isstim(1:end-1)]==0);
Mask_isstim_meanlum_pre3_post3 = (Mask_isstim==0 & [Mask_isstim(2:end); 1]==0 & [Mask_isstim(3:end); ones(2,1)]==0 );

%tmp1_onset = find(tmp1_onset0 & Mask_isstim_meanlum_pre3_post3);

tmp1_onset = ind_run_onset; %in this case, ignore stimuli..

Nonsets = length(tmp1_onset);

WIN = [-2 8];
Fs_USE = Fs1;
t_WIN = [WIN(1)*Fs_USE:WIN(2)*Fs_USE] + WIN(1)*Fs_USE;
Nt = length(t_WIN);
t_WIN2 = t_WIN./Fs;

sz = size(stack4D_ratio);
Run_mat = zeros(Nonsets,sz(1),sz(2),Nt);
Run_mat_dF = zeros(Nonsets,sz(1),sz(2),Nt);
Run_mat_wheell = zeros(Nonsets,Nt);
Run_mat_stim = zeros(Nonsets,Nt);
for count = 1:Nonsets
   ind = tmp1_onset(count) + t_WIN;

   if min(ind)>0 & max(ind)<sz(3)
   Run_mat_wheell(count,:) = tmp1_wheel(ind);
   Run_mat_stim(count,:) = Mask_isstim(ind);
       tmp0 =  double(stack4D_ratio(:,:,ind));
       tmp = tmp0;
       tmp_dF = tmp0;
       ind2 = find(t_WIN2<0 & t_WIN2>-2);
       tmp2 = squeeze(mean(tmp0(:,:,ind2),3));
       for count2  = 1:length(ind)
           tmp(:,:,count2) = (squeeze(tmp0(:,:,count2)) - tmp2)./tmp2;
           tmp_dF(:,:,count2) = (squeeze(tmp0(:,:,count2)) - tmp2);
       end
       Run_mat(count,:,:,:) =tmp;
       Run_mat_dF(count,:,:,:) =tmp_dF;
   end
end

Run_mat_MEAN = squeeze(mean(Run_mat,1));
writetiff(Run_mat_MEAN,[PWD_PRINT2,tmp2_filestr,'Run_mat_onsets_MEANIMJ_NEW.tif']);

Run_mat_dF_MEAN = squeeze(mean(Run_mat_dF,1));
writetiff(Run_mat_dF_MEAN,[PWD_PRINT2,tmp2_filestr,'Run_mat_dF_onsets_MEANIMJ_NEW.tif']);

>>
