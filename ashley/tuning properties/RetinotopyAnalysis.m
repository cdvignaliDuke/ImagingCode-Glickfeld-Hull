SubNum = '516';
date = '141003';
time = '1115';
ImgFolder = '001';
mouse = '516';
fName = '001_000_000';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

% Set current directory to temporary folder on Nuke - cannot analyze data from crash
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder ' - Retinotopy'];
cd(CD);
data_reg = readtiff('Retinotopy_V1.tif');
%%

orig_rate = 30;
final_rate = 3;
down = orig_rate./final_rate;
nON = 100./down;
nOFF = 100./down;
nStim = 6;

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
% max_dF = max(dFoverF_data,[],3);
max_dF = max(dF_data,[],3);
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

%average each trial response per cell
nCells = size(data_TC,2);
dFoverF_meanONbyCell = zeros(nTrials,nCells);
for icell = 1:nCells
    for itrial = 1:nTrials
        tri = data_TC((nON_ind_firsts(itrial):(nON_ind_firsts(itrial)+nON)-1),icell);
        dFoverF_meanONbyCell(itrial,icell) =  mean(tri,1);
    end
end

% resp trace for each trial 
dFoverF_RspbyCellbyTrial = zeros(nON+nOFF,nCells,nTrials);
for icell = 1:nCells
    for itrial = 1:nTrials
        tri = data_TC((nOFF_ind_firsts(itrial):(nON_ind_firsts(itrial)+nON)-1),icell);
        dFoverF_RspbyCellbyTrial(:,icell,itrial) = tri;
    end
end

%remove negative values
minAdd = abs(min(min(dFoverF_RspbyCellbyTrial,[],1)));
dFoverF_RspbyCellbyTrial = bsxfun(@plus, dFoverF_RspbyCellbyTrial,minAdd);
clear minAdd

%% 

% resp matrix for each position
retResp_mat = zeros(nON+nOFF,nCells,nStim,nRep);

for icell = 1:nCells
    start = 1;
    for istim = 1:nStim
        for irep = 1:nRep
            ind = (nStim*(irep-1)+1)+(istim-1);
            retResp_mat(:,icell,istim,irep) = dFoverF_RspbyCellbyTrial(:,icell,ind);
        end
    end
end

%plot
figure;
start = 1;
for icell = 1:5
    for istim = 1:nStim
        pos = istim+(nStim*(start-1));
        subplot(5,6,pos)
        for irep = 1:nRep
            plot(retResp_mat(:,icell,istim,irep));
            hold on
        end
    end
    start = start+1;
end

figure;
col_mat = strvcat('k', 'b', 'c', 'r', 'm', 'g');
for icell = 1:nCells
    subplot(5,6,icell);
    for istim = 1:nStim
        plot(retResp_mat(:,icell,istim,1),col_mat(istim,:));
        hold on
    end
end

