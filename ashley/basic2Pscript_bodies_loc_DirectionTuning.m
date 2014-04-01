%% load MWorks file

mworks = 'data-i023-140317-1616';
load (mworks);

%% Parameters
frame_rate = input.frameImagingRateMs;
orig_rate = 30;
final_rate = 3;
down = orig_rate./final_rate;
nON = 150./down;
nOFF = 150./down;
nStim = 8;
Az = [-30];
El = [15];
Dir = [0 45 90 135 180 225 270 315];
trial_Dir = cell2mat(input.tGratingDirectionDeg);

date = '140317';
mouse = 'G023';

%% IF CONCATINATING DATA SETS USE LOC SECTION OF 'basic2Pscript_ConcatinateFiles.m' 
                            %FROM THIS POINT UNTIL END OF 'Set threshold
                            %for running'
%% Sum-up locomotion data for each trial and calculate speed (Stim ON only)

%create a matrix including all counters
nCounter = floor(((input.nScansOn./frame_rate).*1000)./double(input.speedTimerIntervalMs));
trial_loc_mat = zeros(nCounter,input.trialSinceReset);

for icount = 1:nCounter
    trial_loc_mat(icount, :) = cell2mat_padded(eval(['input.spCounter' num2str(icount)]));
end;

%sum locomotion for each trial (ON stimuli only)
trial_totalloc = sum(trial_loc_mat);

% Calculate speed for each trial in pulses/second
trialon_timeS = input.nScansOn./frame_rate;
trialon_avgspeed_mat = trial_totalloc./trialon_timeS;

%% Set threshold for locomotion and sort trials as running and not running
trialNumbers = 1:length(trialon_avgspeed_mat);
trialon_run = (trialon_avgspeed_mat>4);
runTrials = find(trialon_run);
trialon_LocMatrix = [trialNumbers', trialon_avgspeed_mat',trialon_run'];
trialon_norun = (trialon_avgspeed_mat<4);
norunTrials = find(trialon_norun);
trialon_whichrun = trialon_avgspeed_mat(trialon_run);

%% load 2P imaging data - IF CONCATINATING DATA SETS USE 2P SECTION OF 'basic2Pscript_ConcatinateFiles.m' 
                            %FROM THIS POINT UNTIL END OF 'inde stim types' 
    
data = readrawfile;

%% reshape data
%average signals in time

data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down

%% register
data_avg = mean(data_sub(:,:,300:310),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub
%%
nRep = size(data_reg,3)./((nON+nOFF)*nStim);
%% create dF/F stack

%find off and on frames
nOFF_ind = zeros(1,(nOFF*nStim*nRep));
start = 1;
for iStim = 1:(nRep*nStim)
    nOFF_ind(1, start:start+nOFF-1) = 1+((iStim-1)*(nOFF+nON)):nOFF + ((iStim-1)*(nOFF+nON));
    start = start+nOFF;
end

nON_ind = setdiff(1:size(data_reg,3),nOFF_ind);
nON_avg = mean(data_reg(:,:,nON_ind),3);
nOFF_avg = mean(data_reg(:,:,nOFF_ind),3);

% find on indices for the first frame of each stimulus start period and iti (Off) period
for itrial = 1:(nStim*nRep)
    nON_ind_firsts(itrial) = nON_ind(1+(nON*(itrial-1)));
end
for itrial = 1:(nStim*nRep)
    nOFF_ind_firsts(itrial) = nOFF_ind(1+(nOFF*(itrial-1)));
end

%average nON and nOFF for specific trials
siz = size(data_reg);
trial_nOFF_avg = zeros(siz(1), siz(2) ,size(nOFF_ind_firsts,2));
trial_non_avg = zeros(siz(1), siz(2) ,size(nON_ind_firsts,2));
for itrial = 1:(nStim*nRep)
    trial_nOFF_avg(:,:,itrial) = mean(data_reg(:,:,(nOFF_ind_firsts(itrial): (nOFF_ind_firsts(itrial)+nOFF)-1)),3);
end
for itrial = 1:(nStim*nRep)
    trial_nON_avg(:,:,itrial) = mean(data_reg(:,:,(nON_ind_firsts(itrial): (nON_ind_firsts(itrial)+nON)-1)),3);
end

%dF/F
dF_data = bsxfun(@minus,data_reg, nOFF_avg);
dFoverF_data = bsxfun(@rdivide, dF_data, nOFF_avg);
max_dF = max(dFoverF_data,[],3);
figure; imagesq(max_dF); colormap(gray)

%index stim types
stim_ind = zeros(nStim,nRep);
for idir = 1:nStim
    stim_ind(idir,:) = find(trial_Dir==Dir(idir));
end
%% calculate dF/F on a trial-by-trial basis
trial_dF_data = zeros(siz(1), siz(2) ,(nStim*nRep));
for itrial =  1:(nStim*nRep)
    trial_dF_data(:,:,itrial) = trial_nON_avg(:,:,itrial)-trial_nOFF_avg(:,:,itrial);
end
trial_dFoverF_data = zeros(siz(1), siz(2) ,(nStim*nRep));
for itrial =  1:(nStim*nRep)
    trial_dFoverF_data(:,:,itrial) = trial_dF_data(:,:,itrial)./trial_nOFF_avg(:,:,itrial);
end

%calculate dF/F a stimulus basis and max dF across all stim

stim_ind = zeros(nStim,nRep);
for idir = 1:nStim
    stim_ind(idir,:) = find(trial_Dir==Dir(idir));
end
stim_dFoverF_data = zeros(siz(1),siz(2),nStim);
for idir = 1:nStim
    stim_dFoverF_data(:,:,idir) = mean(trial_dFoverF_data(:,:,stim_ind(idir,:)),3);
end

allstim_max_dF = max(stim_dFoverF_data,[],3);
figure;imagesq(allstim_max_dF); colormap(gray)
%% use max dF/F to find ROIS

bwout = imCellEditInteractive(allstim_max_dF);
mask_cell = bwlabel(bwout);

%timecourses
data_TC = stackGetTimeCourses(dFoverF_data,mask_cell);
figure; tcOffsetPlot(data_TC)
%% average dF/F for each cell per stimulus period
nCell = size(data_TC,2);
for icell = 1:nCell
    for itrial = 1:(nStim*nRep)
        trialon_dFoverFcell(itrial,icell) = mean(data_TC(nON_ind_firsts(itrial):((nON_ind_firsts(itrial)+nON)-1),icell));
    end
end;

%% combine locomotion matrix with frame indices and dF/F means for each cell found

 trialon_LocMatrix = [trialNumbers', trialon_avgspeed_mat',trialon_run',nON_ind_firsts', trial_Dir',trialon_dFoverFcell];
%trialon_LocMatrix = [trialNumbers', trialon_avgspeed_matAll',trialon_run',nON_ind_firsts', trial_Dir_All',trialon_dFoverFcell];

%% which cells are modulated by which stimulus
trialon_dFoverFcell_Dir = [trial_Dir',trialon_dFoverFcell(:,:)];
trialon_dFoverFcell_sort = sortrows(trialon_dFoverFcell_Dir,1);
trialon_dFoverFcell_sorted = trialon_dFoverFcell_sort(:,2:end);
clear trialon_dFoverFcell_sort;

trialon_dFoverFcell_Diravg = zeros(nStim,nCell);

for icell = 1:nCell
    start = 1;
    for istim = 1:nStim
        trialon_dFoverFcell_Diravg(istim,icell) = mean(trialon_dFoverFcell_sorted(start:((start+nRep)-1),icell));
        start = start+nRep;
    end
end

% find cells' preferred orientation
[cell_prefstim_dF_run, cell_prefstim_ind_run] = max(trialon_dFoverFcell_Diravg);

cell_prefstim_type_run = zeros(size(cell_prefstim_ind_run));
cell_prefstim_type_run = Dir(cell_prefstim_ind_run);
cell_pref135_ind = find(cell_prefstim_type_run==135);
cell_pref45_ind = find(cell_prefstim_type_run==45);

%% plot cells dF/F response by orientation 
figure;
for icell = 1:nCell
    plot(Dir, trialon_dFoverFcell_Diravg(:,icell));
    hold on
end

%% plot cells response to locomotion for preferred stimulus

trialon_norun = find(trialon_LocMatrix(:,3)==0);
trialon_run = find(trialon_LocMatrix(:,3)==1);

trialon_LocMatrix_norun = trialon_LocMatrix(trialon_norun,:);
for icell = 1:nCell
    for istim = 1:nStim
        trialon_LocMatrix_norun_ind = find(trialon_LocMatrix_norun(:,5)==Dir(istim));
        if isempty(trialon_LocMatrix_norun_ind) == 0
            trialon_norun_dFoverFcell_Diravg(istim,icell) = mean(trialon_LocMatrix_norun(trialon_LocMatrix_norun_ind,(icell+5)));
        elseif isempty(trialon_LocMatrix_norun_ind) == 1
            trialon_norun_dFoverFcell_Diravg(istim,icell) = 0;
        end
    end
end

trialon_LocMatrix_run = trialon_LocMatrix(trialon_run,:);
for icell = 1:nCell
    for istim = 1:nStim
        trialon_LocMatrix_run_ind = find(trialon_LocMatrix_run(:,5)==Dir(istim));
        if isempty(trialon_LocMatrix_run_ind) == 0
            trialon_run_dFoverFcell_Diravg(istim,icell) = mean(trialon_LocMatrix_run(trialon_LocMatrix_run_ind,(icell+5)));
        elseif isempty(trialon_LocMatrix_run_ind) == 1
            trialon_run_dFoverFcell_Diravg(istim,icell) = 0;
        end
    end
end

% plot run vs. no run trials - all cells
figure;
for icell = 1:nCell
    plot(Dir, trialon_norun_dFoverFcell_Diravg(:,icell),'r');
    hold on
end
hold on
for icell = 1:nCell
    plot(Dir, trialon_run_dFoverFcell_Diravg(:,icell));
    hold on
end
    

%% cells' preferred stim
% find cells' preferred orientation - run
[cell_prefstim_dF_run, cell_prefstim_ind_run] = max(trialon_run_dFoverFcell_Diravg);

cell_prefstim_type_run = Dir(cell_prefstim_ind_run);

% total number of cells that prefer each stimulus condition - run
for istim = 1:nStim
    cell_prefX_ind_run = find(cell_prefstim_type_run==Dir(istim));
    cell_totalpref_run(1,istim) = size(cell_prefX_ind_run,2);
end

% plot run vs. no run trials - cells averaged by preferred stimulus - run
cell_prefstim_avgR_run = zeros(nStim,nStim);
for istimR = 1:nStim
    for istim = 1:nStim
        cell_prefstim_X = find(cell_prefstim_ind_run==istimR);
        cell_prefstim_avgR_run(istim,istimR) = mean(trialon_run_dFoverFcell_Diravg(istim,cell_prefstim_X));
    end
end

% find cells' preferred orientation - no run
[cell_prefstim_dF_norun, cell_prefstim_ind_norun] = max(trialon_norun_dFoverFcell_Diravg);

cell_prefstim_type_norun = Dir(cell_prefstim_ind_norun);

% total number of cells that prefer each stimulus condition - norun
for istim = 1:nStim
    cell_prefX_ind_norun = find(cell_prefstim_type_norun==Dir(istim));
    cell_totalpref_norun(1,istim) = size(cell_prefX_ind_norun,2);
end

% plot run vs. no run trials - cells averaged by preferred stimulus - norun
cell_prefstim_avgR_norun = zeros(nStim,nStim);
for istimR = 1:nStim
    for istim = 1:nStim
        cell_prefstim_X = find(cell_prefstim_ind_norun==istimR);
        cell_prefstim_avgR_norun(istim,istimR) = mean(trialon_norun_dFoverFcell_Diravg(istim,cell_prefstim_X));
    end
end

figure;
for istim = 1:nStim
    plot(Dir, cell_prefstim_avgR_run(istim,:));
    hold on
end
hold on
for istim = 1:nStim
    plot(Dir, cell_prefstim_avgR_norun(istim,:),'r');
    hold on
end

X = find(cell_prefstim_type_run==135);
Y = trialon_norun_dFoverFcell_Diravg(:,X);
nPref = size(Y,2);
figure;
plot(Dir,Y,'-r')
hold on
for icell = 1:nPref
x(:,icell) = find(Y(:,icell)~=0);
plot(Dir,trialon_norun_dFoverFcell_Diravg(x,X),'or')
hold on
plot(Dir,trialon_run_dFoverFcell_Diravg(:,X),'b')

%% save some figs
    sName = strcat('C:\Users\ashley\Documents\mworks\Figs\fig-', num2str(mouse),'-', num2str(date),'-','details', '.pdf');
    epParams = { gcf, sName, ...
             'FileFormat', 'pdf', ...
             'Size', [12 12], ...
             'PrintUI', false };
         exportfig_print(epParams{:});
         close
%% analyze by stimulus type
%find indices for each stim type
stim_mat = zeros(nStim,nRep,(nON+nOFF));
start = 1;
for iRep = 1:nRep
    for iStim = 1:nStim   
        stim_mat(iStim,iRep,:) = 1+((iStim-1)*(nON+nOFF))+((iRep-1)*((nON+nOFF)*nStim)): nON + nOFF + ((iStim-1)*(nON+nOFF))+((iRep-1)*((nON+nOFF)*nStim));
    end
    start= start+nON+nOFF;
end

%plot data
for iCell = 3;
    figure;
    for iStim = 1:nStim
        subplot(2,3,iStim)
        rep_mat = zeros(nON+nOFF,nRep);
        for iRep = 1:nRep
            plot(1-nOFF:nON,data_TC(squeeze(stim_mat(iStim,iRep,:))',iCell), 'k');
            hold on
            rep_mat(:,iRep) = data_TC(squeeze(stim_mat(iStim,iRep,:))',iCell);
        end
        plot(1-nOFF:nON, mean(rep_mat,2), 'r');
        hold on
        ylim([min(data_TC(:,iCell),[],1) max(data_TC(:,iCell),[],1)])
        stim_Az = rem(iStim,size(Az,2));
        if stim_Az == 0
            stim_Az = size(Az,2);
        end
        stim_El= ceil(iStim./size(Az,2));
        title(['Az = ' num2str(Az(1,stim_Az)) ' El = ' num2str(El(1,stim_El))]);
    end
end
