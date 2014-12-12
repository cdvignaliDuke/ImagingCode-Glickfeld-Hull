%% Parameters
% frame_rate = input.frameImagingRateMs;
orig_rate = 30;
final_rate = 3;
down = orig_rate./final_rate;
nON = 100./down;
nOFF = 100./down;
nStim = 6;
Az = [0 15 30];
El = [10 -10];
position = 1:6;
FSpos = 3;

%% reshape data
%average signals in time

data_down = stackGroupProject(data,down);
clear data

%remove negative data (by addition)
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down

% register
data_avg = mean(data_sub(:,:,300:310),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

%save data_reg
writetiff(data_reg, 'Retinotopy_V1');

%registered image
data_avg = mean(data_reg(:,:,:),3);
figure; imagesq(data_avg); colormap(gray)

%%
nRep = size(data_reg,3)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);
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


%% which fields of view or single cells are modulated by which stimulus?
trialAZ = (cell2mat(input.tGratingAzimuthDeg))';
trialEL = (cell2mat(input.tGratingElevationDeg))';
pos_mat = [trialAZ trialEL];

%write a for loop to name and fill these variables...
%Index for where each type (position) of stimulus is by trial. 
pos1_ind = find((pos_mat(:,1)== Az(1)) & (pos_mat(:,2) == El(1)));
pos2_ind = find((pos_mat(:,1)== Az(2)) & (pos_mat(:,2) == El(1)));
pos3_ind = find((pos_mat(:,1)== Az(3)) & (pos_mat(:,2) == El(1)));
pos4_ind = find((pos_mat(:,1)== Az(1)) & (pos_mat(:,2) == El(2)));
pos5_ind = find((pos_mat(:,1)== Az(2)) & (pos_mat(:,2) == El(2)));
pos6_ind = find((pos_mat(:,1)== Az(3)) & (pos_mat(:,2) == El(2)));
% pos7_ind = find((pos_mat(:,1)== Az(2)) & (pos_mat(:,2) == El(2)));
% pos8_ind = find((pos_mat(:,1)== Az(3)) & (pos_mat(:,2) == El(2)));
% pos9_ind = find((pos_mat(:,1)== Az(4)) & (pos_mat(:,2) == El(2)));
% pos10_ind = find((pos_mat(:,1)== Az(5)) & (pos_mat(:,2) == El(2)));

%find average response to a particular stimulus - full field
pos1_respavg = mean(trialon_dFoverFall(:,:,pos1_ind),3);
pos2_respavg = mean(trialon_dFoverFall(:,:,pos2_ind),3);
pos3_respavg = mean(trialon_dFoverFall(:,:,pos3_ind),3);
pos4_respavg = mean(trialon_dFoverFall(:,:,pos4_ind),3);
pos5_respavg = mean(trialon_dFoverFall(:,:,pos5_ind),3);
pos6_respavg = mean(trialon_dFoverFall(:,:,pos6_ind),3);
% pos7_respavg = mean(trialon_dFoverFall(:,:,pos7_ind),3);
% pos8_respavg = mean(trialon_dFoverFall(:,:,pos8_ind),3);
% pos9_respavg = mean(trialon_dFoverFall(:,:,pos9_ind),3);
% pos10_respavg = mean(trialon_dFoverFall(:,:,pos10_ind),3);

%find max full-field response to stimulus
figure;imagesq(pos1_respavg);colormap gray
imagesq(pos2_respavg);colormap gray
imagesq(pos3_respavg);colormap gray
imagesq(pos4_respavg);colormap gray
imagesq(pos5_respavg);colormap gray
imagesq(pos6_respavg);colormap gray

Retinotopy = figure; 
subplot(2, 3, 1);
imagesq(pos1_respavg); title('Az:-15 El:15'); colormap gray; caxis([-.5 .5]); colorbar;
hold on
subplot(2,3,2);
imagesq(pos2_respavg); title('Az:0 El:15'); colormap gray; caxis([-.5 .5]); colorbar;
hold on 
subplot(2,3,3);
imagesq(pos3_respavg); title('Az:15 El:15'); colormap gray; caxis([-.5 .5]); colorbar;
subplot(2,3,4);
imagesq(pos4_respavg); title('Az:-15 El:0'); colormap gray; caxis([-.5 .5]); colorbar;
subplot(2,3,5);
imagesq(pos5_respavg); title('Az:0 El:0'); colormap gray; caxis([-.5 .5]); colorbar;
subplot(2,3,6);
imagesq(pos6_respavg); title('Az:15 El:0'); colormap gray; caxis([-.5 .5]); colorbar;



%save fig
saveas(Retinotopy,'Retinotopy.fig')

%% average dF/F for each cell per stimulus period
nCell = size(data_TC,2);
for icell = 1:nCell
    for itrial = 1:nTrials
        trialon_dFoverFcell(itrial,icell) = mean(data_TC(nON_ind_firsts(itrial):((nON_ind_firsts(itrial)+nON)-1),icell));
    end
end;
%find average response to a particular stimulus - by cell
pos1_cellavg = zeros(1,nCell);
for icell = 1:nCell
    pos1_cellavg(:,icell) = mean(trialon_dFoverFcell(pos1_ind,icell));
end
pos2_cellavg = zeros(1,nCell);
for icell = 1:nCell
    pos2_cellavg(:,icell) = mean(trialon_dFoverFcell(pos2_ind,icell));
end
pos3_cellavg = zeros(1,nCell);
for icell = 1:nCell
    pos3_cellavg(:,icell) = mean(trialon_dFoverFcell(pos3_ind,icell));
end
pos4_cellavg = zeros(1,nCell);
for icell = 1:nCell
    pos4_cellavg(:,icell) = mean(trialon_dFoverFcell(pos4_ind,icell));
end
pos5_cellavg = zeros(1,nCell);
for icell = 1:nCell
    pos5_cellavg(:,icell) = mean(trialon_dFoverFcell(pos5_ind,icell));
end
pos6_cellavg = zeros(1,nCell);
for icell = 1:nCell
    pos6_cellavg(:,icell) = mean(trialon_dFoverFcell(pos6_ind,icell));
end
%matrix for all trials by all cells
 pos_cellavg_mat = [pos1_cellavg; pos2_cellavg; pos3_cellavg; pos4_cellavg; pos5_cellavg; pos6_cellavg];
 
%strongest fluorescence by trial for each cell
pos_max_ind = zeros(1,nCell);
for icell = 1:nCell
    pos_max_ind(1,icell) = find(pos_cellavg_mat(:,icell) == max(pos_cellavg_mat(:,icell)));
end

pos_pref_totalcells = zeros(1,nStim);
for ipos = 1:nStim
    pos_pref_ind = find(pos_max_ind == ipos);
    pos_pref_totalcells(:,ipos) = size(pos_pref_ind,2);    
end

%index for cells that prefer the posititon used in the Flashing Stim experiment
FSpos_ind = find(pos_max_ind == FSpos);

%save variables
save('analysis')
