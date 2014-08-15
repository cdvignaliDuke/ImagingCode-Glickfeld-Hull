%% load data
SubNum = '001';
date = '140815';
time = '006';
ImgFolder = '006';
mouse = 'AW01';

% load MWorks file
CD = ['Z:\2P imaging\MWorks\' mouse '\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

% Set current directory to temporary folder on Nuke - cannot analyze data from crash
CD = ['D:\Ashley_temp' '\' ImgFolder];
cd(CD);

%% Parameters
% frame_rate = input.frameImagingRateMs;
orig_rate = 30;
final_rate = 3;
down = orig_rate./final_rate;
nON = 150./down;
nOFF = 150./down;
nStim = 1;
Az = [30];
El = [10];

%% load 2P imaging data 

%set the number of parallel streams (I think we have up to 8 cores)
tic
ncores = 5;
parpool = ncores;
%set the size of the batch that each stream will read
nbatch = 600;
%create a parallel loop to load the data
fName = '006_000_000';
parfor i=1:5
    data(:,:,:,:,i) =sbxread(fName,(i-1)*nbatch,nbatch);
end
%reshape z to make a 3d stack
pmt = 1; %1 = green 2 = red
data = squeeze(data(pmt,:,:,:,:));
siz = size(data);
reshape(data,siz(1),siz(2),siz(3)*ncores);toc

    
%data = readrawfile;

%% reshape data
%average signals in time

data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down

% register
data_avg = mean(data_sub(:,:,20:30),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

%registered image
data_avg = mean(data_reg(:,:,:),3);
figure; imagesq(data_avg); colormap(gray)

%%chunk data set (1000 frames apiece) to load all at once without taking up
    %too much memory
data_reg1 = data_reg;

%chunk 2
NumFrames2 = 2000;
data = sbxread(fName,(NumFrames+1),NumFrames);
% data will be a 4D matrix [#pmts rows(y) columns(x) NumFrames]
%choose pmt to look at 
pmt = 1; %1 = green 2 = red
data = data(pmt,:,:,:);
data = squeeze(data);
data_down = stackGroupProject(data,down);
clear data
%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down
% register
data_avg = mean(data_sub(:,:,20:30),3);
[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub
data_reg2 = data_reg;
% 
% %chunk 3
% NumFrames3 = 3000;
% data = sbxread(fName,(NumFrames2+1),NumFrames);
% % data will be a 4D matrix [#pmts rows(y) columns(x) NumFrames]
% %choose pmt to look at 
% pmt = 1; %1 = green 2 = red
% data = data(pmt,:,:,:);
% data = squeeze(data);
% data_down = stackGroupProject(data,down);
% clear data
% %remove negative data (by addition)
% data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
% clear data_down
% % register
% data_avg = mean(data_sub(:,:,20:30),3);
% [out data_reg] = stackRegister(data_sub, data_avg);
% clear data_sub
% data_reg3 = data_reg;
% 
% %chunk 4
% NumFrames4 = 4000;
% data = sbxread(fName,(NumFrames3+1),NumFrames);
% % data will be a 4D matrix [#pmts rows(y) columns(x) NumFrames]
% %choose pmt to look at 
% pmt = 1; %1 = green 2 = red
% data = data(pmt,:,:,:);
% data = squeeze(data);
% data_down = stackGroupProject(data,down);
% clear data
% %remove negative data (by addition)
% data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
% clear data_down
% % register
% data_avg = mean(data_sub(:,:,20:30),3);
% [out data_reg] = stackRegister(data_sub, data_avg);
% clear data_sub
% data_reg4 = data_reg;
% 
% %chunk 5
% NumFrames5 = 5000;
% data = sbxread(fName,(NumFrames4+1),NumFrames);
% % data will be a 4D matrix [#pmts rows(y) columns(x) NumFrames]
% %choose pmt to look at 
% pmt = 1; %1 = green 2 = red
% data = data(pmt,:,:,:);
% data = squeeze(data);
% data_down = stackGroupProject(data,down);
% clear data
% %remove negative data (by addition)
% data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
% clear data_down
% % register
% data_avg = mean(data_sub(:,:,20:30),3);
% [out data_reg] = stackRegister(data_sub, data_avg);
% clear data_sub
% data_reg5 = data_reg;
% 
% %chunk 6
% NumFrames6 = 6000;
% data = sbxread(fName,(NumFrames5+1),NumFrames);
% % data will be a 4D matrix [#pmts rows(y) columns(x) NumFrames]
% %choose pmt to look at 
% pmt = 1; %1 = green 2 = red
% data = data(pmt,:,:,:);
% data = squeeze(data);
% data_down = stackGroupProject(data,down);
% clear data
% %remove negative data (by addition)
% data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
% clear data_down
% % register
% data_avg = mean(data_sub(:,:,20:30),3);
% [out data_reg] = stackRegister(data_sub, data_avg);
% clear data_sub
% data_reg6 = data_reg;
% 
% %chunk 7
% NumFrames7 = 7000;
% data = sbxread(fName,(NumFrames6+1),NumFrames);
% % data will be a 4D matrix [#pmts rows(y) columns(x) NumFrames]
% %choose pmt to look at 
% pmt = 1; %1 = green 2 = red
% data = data(pmt,:,:,:);
% data = squeeze(data);
% data_down = stackGroupProject(data,down);
% clear data
% %remove negative data (by addition)
% data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
% clear data_down
% % register
% data_avg = mean(data_sub(:,:,20:30),3);
% [out data_reg] = stackRegister(data_sub, data_avg);
% clear data_sub
% data_reg7 = data_reg;
% 
% %chunk 8
% NumFrames8 = 8000;
% data = sbxread(fName,(NumFrames7+1),NumFrames);
% % data will be a 4D matrix [#pmts rows(y) columns(x) NumFrames]
% %choose pmt to look at 
% pmt = 1; %1 = green 2 = red
% data = data(pmt,:,:,:);
% data = squeeze(data);
% data_down = stackGroupProject(data,down);
% clear data
% %remove negative data (by addition)
% data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
% clear data_down
% % register
% data_avg = mean(data_sub(:,:,20:30),3);
% [out data_reg] = stackRegister(data_sub, data_avg);
% clear data_sub
% data_reg8 = data_reg;
% 
% %chunk 9
% NumFrames9 = 9000;
% data = sbxread(fName,(NumFrames8+1),NumFrames);
% % data will be a 4D matrix [#pmts rows(y) columns(x) NumFrames]
% %choose pmt to look at 
% pmt = 1; %1 = green 2 = red
% data = data(pmt,:,:,:);
% data = squeeze(data);
% data_down = stackGroupProject(data,down);
% clear data
% %remove negative data (by addition)
% data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
% clear data_down
% % register
% data_avg = mean(data_sub(:,:,20:30),3);
% [out data_reg] = stackRegister(data_sub, data_avg);
% clear data_sub
% data_reg9 = data_reg;
% 
% %put in a for loop for this chunking...
% 
% %concatinate 
% data_reg = cat(3,data_reg1,data_reg2,data_reg3,data_reg4,data_reg5,data_reg6,data_reg7,data_reg8,data_reg9);
% clear data_reg1 data_reg2 data_reg3 data_reg4 data_reg5 data_reg6 data_reg7 data_reg8 data_reg9

%%
nRep = size(data_reg,3)./((nON+nOFF)*nStim);
nTrials = 1:(nStim.*nRep);
%% create dF/F stack

%find off and on frames
nOFF_ind = zeros(1,(nOFF*nStim*nRep));
start = 1;
for iStim = 1:(nRep*nStim)
    nOFF_ind(1, start:start+nOFF-1) = 1+((iStim-1)*(nOFF+nON)):nOFF + ((iStim-1)*(nOFF+nON));
    start = start+nOFF;
end

nON_ind = setdiff(1:size(data_reg,3),nOFF_ind);
nON_avg = mean(x(:,:,nON_ind),3);
nOFF_avg = mean(x(:,:,nOFF_ind),3);

%dF/F
dF_data = bsxfun(@minus,data_reg, nOFF_avg);
dFoverF_data = bsxfun(@rdivide, dF_data, nOFF_avg);
max_dF = max(dFoverF_data,[],3);
figure; imagesq(max_dF); colormap(gray)

%% use max dF/F to find ROIS

bwout = imCellEditInteractive(max_dF);
mask_cell = bwlabel(bwout);

%timecourses
data_TC = stackGetTimeCourses(dFoverF_data,mask_cell);
figure; tcOffsetPlot(data_TC)
%%
% find on indices for the first frame of each stimulus start period and iti (Off) period
for itrial = 1:(nStim*nRep)
    nON_ind_firsts(itrial) = nON_ind(1+(nON*(itrial-1)));
end
for itrial = 1:(nStim*nRep)
    nOFF_ind_firsts(itrial) = nOFF_ind(1+(nOFF*(itrial-1)));
end

%% average dF/F for each field of view per stimulus period
trialon_dFoverFall = zeros((size(dFoverF_data,1)),(size(dFoverF_data,2)),nTrials);
for itrial = 1:nTrials
    tri = dFoverF_data(:,:,(nON_ind_firsts(itrial):(nON_ind_firsts(itrial)+nON)-1));
    trialon_dFoverFall (:,:,itrial) =  mean(tri,3);
end

%% average dF/F for each cell per stimulus period
nCell = size(data_TC,2);
for icell = 1:nCell
    for itrial = 1:nTrials
        trialon_dFoverFcell(itrial,icell) = mean(data_TC(nON_ind_firsts(itrial):((nON_ind_firsts(itrial)+nON)-1),icell));
    end
end;

%% combine locomotion matrix with frame indices and dF/F means for each cell found

%  trialon_LocMatrix = [trialNumbers', trialon_avgspeed_mat',trialon_run',nON_ind_firsts', trial_Dir',trialon_dFoverFcell];
%trialon_LocMatrix = [trialNumbers', trialon_avgspeed_matAll',trialon_run',nON_ind_firsts', trial_Dir_All',trialon_dFoverFcell];

%% which fields of view or single cells are modulated by which stimulus?
trialAZ = (cell2mat(input.tGratingAzimuthDeg))';
trialEL = (cell2mat(input.tGratingElevationDeg))';
pos_mat = [trialAZ trialEL];

%write a for loop to name and fill these variables...
pos1_ind = find((pos_mat(:,1)== Az(1)) & (pos_mat(:,2) == El(1)));
pos2_ind = find((pos_mat(:,1)== Az(2)) & (pos_mat(:,2) == El(1)));
pos3_ind = find((pos_mat(:,1)== Az(3)) & (pos_mat(:,2) == El(1)));
pos4_ind = find((pos_mat(:,1)== Az(4)) & (pos_mat(:,2) == El(1)));
pos5_ind = find((pos_mat(:,1)== Az(5)) & (pos_mat(:,2) == El(1)));
pos6_ind = find((pos_mat(:,1)== Az(1)) & (pos_mat(:,2) == El(2)));
pos7_ind = find((pos_mat(:,1)== Az(2)) & (pos_mat(:,2) == El(2)));
pos8_ind = find((pos_mat(:,1)== Az(3)) & (pos_mat(:,2) == El(2)));
pos9_ind = find((pos_mat(:,1)== Az(4)) & (pos_mat(:,2) == El(2)));
pos10_ind = find((pos_mat(:,1)== Az(5)) & (pos_mat(:,2) == El(2)));

%find average response to a particular stimulus - full field
pos1_respavg = mean(trialon_dFoverFall(:,:,pos1_ind),3);
pos2_respavg = mean(trialon_dFoverFall(:,:,pos2_ind),3);
pos3_respavg = mean(trialon_dFoverFall(:,:,pos3_ind),3);
pos4_respavg = mean(trialon_dFoverFall(:,:,pos4_ind),3);
pos5_respavg = mean(trialon_dFoverFall(:,:,pos5_ind),3);
pos6_respavg = mean(trialon_dFoverFall(:,:,pos6_ind),3);
pos7_respavg = mean(trialon_dFoverFall(:,:,pos7_ind),3);
pos8_respavg = mean(trialon_dFoverFall(:,:,pos8_ind),3);
pos9_respavg = mean(trialon_dFoverFall(:,:,pos9_ind),3);
pos10_respavg = mean(trialon_dFoverFall(:,:,pos10_ind),3);

%find max full-field response to stimulus
pos1_val = sum(sum(pos1_respavg));
pos2_val = sum(sum(pos2_respavg));
pos3_val = sum(sum(pos3_respavg));
pos4_val = sum(sum(pos4_respavg));
pos5_val = sum(sum(pos5_respavg));
pos6_val = sum(sum(pos6_respavg));
pos7_val = sum(sum(pos7_respavg));
pos8_val = sum(sum(pos8_respavg));
pos9_val = sum(sum(pos9_respavg));
pos10_val = sum(sum(pos10_respavg));
avg_responses_all = [pos1_val pos2_val pos3_val pos4_val pos5_val pos6_val pos7_val pos8_val pos9_val pos10_val]

    %plot that
hist(avg_responses_all);

% trialon_dFoverFcell_Pos = [trial_Dir',trialon_dFoverFcell(:,:)];
% trialon_dFoverFcell_sort = sortrows(trialon_dFoverFcell_Dir,1);
% trialon_dFoverFcell_sorted = trialon_dFoverFcell_sort(:,2:end);
% clear trialon_dFoverFcell_sort;

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
