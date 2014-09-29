%% sort data_TC by visual and auditory+visual trials

Block2ON_005 = double(cell2mat(input.tBlock2TrialNumber));
Block2ON_006 = double(cell2mat(input.tBlock2TrialNumber));
Block2ON = cat(2,Block2ON_005,Block2ON_006);
Block2ON = Block2ON(:,1:119);

% Block2ON = Block2ON(:,1:end-1);

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
A_ind_sub = A_ind(find((A_ind>10 & A_ind<60)|(A_ind>69)));
%A_ind_sub = A_ind(find((A_ind<10)|(60<A_ind<69)));
V_ind = find(Block2ON == 0);
V_ind_sub = V_ind(find((V_ind>10 & V_ind<60)|(V_ind>69)));
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

figure; plot(Cell_3visstim_mean_A(:,visualresp_ind),'r')
hold on; plot(Cell_3visstim_mean_V(:,visualresp_ind), 'g')

figure;
start = 1;
for icell = visualresp_ind
    subplot(3,3,start);
    plot(Cell_3visstim_mean_A(:,icell), 'r');
    hold on;
    plot(Cell_3visstim_mean_V(:,icell), 'g');
    axis([0 42 -0.1 0.1]);
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
ylabel('dF/F');
vline(10,'k-');
vline(20.5,'k-');
vline(31,'k-');
title(['Average All Cells - Auditory vs. Visual']);


%plot mask_cell for these cells
figure;
start = 1;
for icell = visualresp_ind
    subplot(3,4,start);
    mask_cell_new = zeros(size(mask_cell));
    mask_cell_new(find(mask_cell==icell)) = 1;
    imagesc(mask_cell_new);
    title(['Cell ' num2str(icell)]);
    start = start + 1;
end

%% Control - randomized trials.

RandomON = randi([0 1],1,tTrials);

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
axis([0 42 -0.02 0.04]);
ylabel('dF/F');
vline(10,'k-');
vline(20.5,'k-');
vline(31,'k-');
title(['Average All Cells -Trials separated randomly']);



figure;
start = 1;
for icell = 1:6
    subplot(2,3,start);
    plot(Cell_3visstim_mean_A(:,icell), 'm');
    hold on;
    plot(Cell_3visstim_mean_V(:,icell), 'c');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end

figure;
start = 1;
for icell = 7:12
    subplot(2,3,start);
    plot(Cell_3visstim_mean_A(:,icell), 'm');
    hold on;
    plot(Cell_3visstim_mean_V(:,icell), 'c');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end

figure;
start = 1;
for icell = 13:18
    subplot(2,3,start);
    plot(Cell_3visstim_mean_A(:,icell), 'm');
    hold on;
    plot(Cell_3visstim_mean_V(:,icell), 'c');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end

figure;
start = 1;
for icell = 19:24
    subplot(2,3,start);
    plot(Cell_3visstim_mean_A(:,icell), 'm');
    hold on;
    plot(Cell_3visstim_mean_V(:,icell), 'c');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end

figure;
start = 1;
for icell = 25:30
    subplot(2,3,start);
    plot(Cell_3visstim_mean_A(:,icell), 'm');
    hold on;
    plot(Cell_3visstim_mean_V(:,icell), 'c');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end

figure;
start = 1;
for icell = 31:36
    subplot(2,3,start);
    plot(Cell_3visstim_mean_A(:,icell), 'm');
    hold on;
    plot(Cell_3visstim_mean_V(:,icell), 'c');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end

figure;
start = 1;
for icell = 36:nCells
    subplot(2,3,start);
    plot(Cell_3visstim_mean_A(:,icell), 'm');
    hold on;
    plot(Cell_3visstim_mean_V(:,icell), 'c');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end

