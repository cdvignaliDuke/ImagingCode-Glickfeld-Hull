%save path
ImgFolder = '002+003+004';
Save = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\Rsp2VisStim'];
cd(Save)

%% *********behavior experiment***************
%% sort data_TC by success & ignore, visual and auditory 
Block2ONSuccessAndIgnore = Block2ON(:,SuccessANDIgnoreTrials_ind);

nCells = size(data_TC,2);
Cell_preLeverDownplus3visstim_mat = zeros(L,tTrials,nCells);
Cell_preLeverDownplus3visstim_trialmean = zeros(L,nCells);
for icell = 1:nCells
    for itrial = 1:tTrials
        Cell_preLeverDownplus3visstim_mat(:,itrial,icell) = data_TC(1+(42.*(itrial-1)):42.*itrial,icell);
    end
    Cell_preLeverDownplus3visstim_trialmean(:,icell) = mean(Cell_preLeverDownplus3visstim_mat(:,:,icell),2);
end

A_ind = find(Block2ONSuccessAndIgnore == 1);
V_ind = find(Block2ONSuccessAndIgnore == 0);

Cell_3visstim_mat_A = Cell_preLeverDownplus3visstim_mat(:,A_ind,:);
for icell = 1:nCells
    Cell_3visstim_mean_A(:,icell) = mean(Cell_3visstim_mat_A(:,:,icell),2);
end
Cell_3visstim_mat_V = Cell_preLeverDownplus3visstim_mat(:,V_ind,:);
for icell = 1:nCells
    Cell_3visstim_mean_V(:,icell) = mean(Cell_3visstim_mat_V(:,:,icell),2);
end

figure;
start = 1;
for icell = 1:nCells
    subplot(6,5,start);
    plot(Cell_3visstim_mean_A(:,icell), 'r');
    hold on;
    plot(Cell_3visstim_mean_V(:,icell), 'g');
    axis([0 42 -0.1 0.3]);
    vline(10,'k-');
    vline(20.5,'k-');
    vline(31,'k-');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end

AllCells_3visstim_mean_A = mean(Cell_3visstim_mean_A,2);
AllCells_3visstim_mean_V = mean(Cell_3visstim_mean_V,2);

figure; plot(AllCells_3visstim_mean_A,'r')
hold on; plot(AllCells_3visstim_mean_V, 'g')
axis([0 42 -0.02 0.2]);
ylabel('dF/F');
vline(10,'k-');
vline(20.5,'k-');
vline(31,'k-');
title(['Average All Cells - Auditory vs. Visual']);

%% compare first vs. last trial
RspFirstTrial_A = mean(Cell_preLeverDownplus3visstim_mat(:,[1 5 7],visualresp_ind),3);
RspLastTrial_A = mean(Cell_preLeverDownplus3visstim_mat(:,[57 61 64],visualresp_ind),3);

figure; plot(RspFirstTrial_A,'m')
hold on; plot(RspLastTrial_A, 'b')
axis([0 42 -0.02 0.8]);
ylabel('dF/F');
vline(10,'k-');
vline(20.5,'k-');
vline(31,'k-');
title(['Average Responisve Cells - First vs. Last, Auditory']);

RspFirstTrial_V = mean(Cell_preLeverDownplus3visstim_mat(:,[2 3 4],visualresp_ind),3);
RspLastTrial_V = mean(Cell_preLeverDownplus3visstim_mat(:,[62 63 65],visualresp_ind),3);

figure; plot(RspFirstTrial_V,'m')
hold on; plot(RspLastTrial_V, 'b')
axis([0 42 -0.02 0.8]);
ylabel('dF/F');
vline(10,'k-');
vline(20.5,'k-');
vline(31,'k-');
title(['Average Responisve Cells - First vs. Last, Visual']);

%% compare success vs. ignore trials
SuccessTrials_ind_A = find(Block2ON(:,SuccessTrials_ind)==1);
IgnoreTrials_ind_A = find(Block2ON(:,IgnoreTrials_ind)==1);

RspSucTrial_A = mean(Cell_preLeverDownplus3visstim_mat(:,SuccessTrials_ind_A,visualresp_ind),3);
RspIgTrial_A = mean(Cell_preLeverDownplus3visstim_mat(:,IgnoreTrials_ind_A,visualresp_ind),3);

figure; plot(RspSucTrial_A, 'k');
hold on; plot(RspIgTrial_A, 'm');
axis([0 42 -0.02 0.8]);
ylabel('dF/F');
vline(10,'k-');
vline(20.5,'k-');
vline(31,'k-');
title(['Average Responisve Cells - Success vs Ignore, Auditory']);

SuccessTrials_ind_V = find(Block2ON(:,SuccessTrials_ind)==0);
IgnoreTrials_ind_V = find(Block2ON(:,IgnoreTrials_ind)==0);

RspSucTrial_V = mean(Cell_preLeverDownplus3visstim_mat(:,SuccessTrials_ind_V,visualresp_ind),3);
RspIgTrial_V = mean(Cell_preLeverDownplus3visstim_mat(:,IgnoreTrials_ind_V,visualresp_ind),3);

figure; plot(RspSucTrial_V, 'k');
hold on; plot(RspIgTrial_V, 'm');
axis([0 42 -0.02 0.8]);
ylabel('dF/F');
vline(10,'k-');
vline(20.5,'k-');
vline(31,'k-');
title(['Average Responisve Cells - Success vs Ignore, Visual']);


AllCells_AvgSuccessTrials_A = mean(RspSucTrial_A,2);
AllCells_AvgIgnoreTrials_A = mean(RspIgTrial_A,2);
AllCells_AvgSuccessTrials_V = mean(RspSucTrial_V,2);
AllCells_AvgIgnoreTrials_V = mean(RspIgTrial_V,2);

figure; plot(AllCells_AvgSuccessTrials_A, 'b')
hold on; plot(AllCells_AvgIgnoreTrials_A, 'y')
hold on; plot(AllCells_AvgSuccessTrials_V, 'k')
hold on; plot(AllCells_AvgIgnoreTrials_V, 'm')
axis([0 42 -0.02 0.1]);
ylabel('dF/F');
vline(10,'k-');
vline(20.5,'k-');
vline(31,'k-');
title(['Average all cells, all trials, Success vs Ignore,Visual and auditory']);


%% ***********Fake mouse experiment*************** 
%% sort data_TC by visual and auditory+visual trials

Block2ON_008 = cell2mat(input.tBlock2TrialNumber);
Block2ON_008 = Block2ON_008(:,1:end-1);
Block2ON_009 = cell2mat(input.tBlock2TrialNumber);
Block2ON_009 = Block2ON_009(:,1:end-1);
Block2ON = cat(2,Block2ON_008,Block2ON_009);


nCells = size(data_TC,2);
Cell_preLeverDownplus3visstim_mat = zeros(L,tTrials,nCells);
Cell_preLeverDownplus3visstim_trialmean = zeros(L,nCells);
for icell = 1:nCells
    for itrial = 1:tTrials
        Cell_preLeverDownplus3visstim_mat(:,itrial,icell) = data_TC(1+(42.*(itrial-1)):42.*itrial,icell);
    end
    Cell_preLeverDownplus3visstim_trialmean(:,icell) = mean(Cell_preLeverDownplus3visstim_mat(:,:,icell),2);
end

A_ind = find(Block2ON == 1);
% A_ind_sub = find(Block2ON_test == 1);
%A_ind_sub = A_ind(find((A_ind<10)|(60<A_ind<69)));
V_ind = find(Block2ON == 0);
% V_ind_sub = find(Block2ON_test == 0); 
%V_ind_sub = V_ind(find((V_ind<10)|(60<V_ind>69)));

Cell_3visstim_mat_A = Cell_preLeverDownplus3visstim_mat(:,A_ind,:);
for icell = 1:nCells
    Cell_3visstim_mean_A(:,icell) = mean(Cell_3visstim_mat_A(:,:,icell),2);
end
Cell_3visstim_mat_V = Cell_preLeverDownplus3visstim_mat(:,V_ind,:);
for icell = 1:nCells
    Cell_3visstim_mean_V(:,icell) = mean(Cell_3visstim_mat_V(:,:,icell),2);
end

% Cell_3visstim_mat_A = Cell_preLeverDownplus3visstim_mat(:,A_ind_sub,:);
% for icell = 1:nCells
%     Cell_3visstim_mean_A(:,icell) = mean(Cell_3visstim_mat_A(:,:,icell),2);
% end
% Cell_3visstim_mat_V = Cell_preLeverDownplus3visstim_mat(:,V_ind_sub,:);
% for icell = 1:nCells
%     Cell_3visstim_mean_V(:,icell) = mean(Cell_3visstim_mat_V(:,:,icell),2);
% end

figure; plot(Cell_3visstim_mean_A,'r')
hold on; plot(Cell_3visstim_mean_V, 'g')

figure;
start = 1;
for icell = 1:nCells
    subplot(5,5,start);
    plot(Cell_3visstim_mean_A(:,icell), 'r');
    hold on;
    plot(Cell_3visstim_mean_V(:,icell), 'g');
    axis([0 42 -0.1 0.2]);
    vline(10,'k-');
    vline(20.5,'k-');
    vline(31,'k-');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end

AllCells_3visstim_mean_A = mean(Cell_3visstim_mean_A,2);
AllCells_3visstim_mean_V = mean(Cell_3visstim_mean_V,2);

figure; plot(AllCells_3visstim_mean_A,'r')
hold on; plot(AllCells_3visstim_mean_V, 'g')
axis([0 42 -0.02 0.04]);
vline(10,'k-');
vline(20.5,'k-');
vline(31,'k-');
title(['Average All Cells - Auditory vs. Visual']);


%plot mask_cell for these cells
figure;
start = 1;
for icell = 19
    subplot(3,4,start);
    mask_cell_new = zeros(size(mask_cell));
    mask_cell_new(find(mask_cell==icell)) = 1;
    imagesc(mask_cell_new);
    title(['Cell ' num2str(icell)]);
    start = start + 1;
end

%% Control - randomized trials.

% %fake mouse data
% RandomON = randi([0 1],1,tTrials);

%behavior data
tTrialsSuccessAndIgnore = size(SuccessANDIgnoreTrials_ind,2);
RandomON = randi([0 1],1,tTrialsSuccessAndIgnore);

R1_ind = find(RandomON == 1);
R0_ind = find(RandomON == 0);
Cell_3visstim_mat_R1 = Cell_preLeverDownplus3visstim_mat(:,R1_ind,:);
for icell = visualresp_ind 
    Cell_3visstim_mean_R1(:,icell) = mean(Cell_3visstim_mat_R1(:,:,icell),2);
end
Cell_3visstim_mat_R0 = Cell_preLeverDownplus3visstim_mat(:,R0_ind,:);
for icell = visualresp_ind
    Cell_3visstim_mean_R0(:,icell) = mean(Cell_3visstim_mat_R0(:,:,icell),2);
end


figure; plot(Cell_3visstim_mean_R1,'m')
hold on; plot(Cell_3visstim_mean_R0, 'c')

AllCells_3visstim_mean_R1 = mean(Cell_3visstim_mean_R1,2);
AllCells_3visstim_mean_R0 = mean(Cell_3visstim_mean_R0,2);

figure; plot(AllCells_3visstim_mean_R1,'m')
hold on; plot(AllCells_3visstim_mean_R0, 'c')
axis([0 42 -0.02 0.3]);
ylabel('dF/F');
vline(10,'k-');
vline(20.5,'k-');
vline(31,'k-');
title(['Average All Cells -Trials separated randomly']);



