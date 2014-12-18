%load data
SubNum = '607';
date = '141209';
time = '1807';
ImgFolder = '004';
mouse = 'AW07';
fName = '004_000_000';

% load MWorks file
% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

% Set current directory to temporary folder on Nuke - cannot analyze data from crash
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
data = readtiff('DirectionTuning_V1.tif');

%% Parameters
orig_rate = 30;
final_rate = 3;
down = orig_rate./final_rate;
nON = (input.nScansOn)./down;
nOFF = (input.nScansOff)./down;
nStim = input.gratingDirectionStepN;

%%

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
writetiff(data_reg, 'DirectionTuning_V1');

%%
nRep = size(data_reg,3)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);
% DirectionDeg = cell2mat(input.tGratingDirectionDeg);
% Dirs = unique(DirectionDeg);
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
figure; imagesq(max_dF); colormap(gray)

%% use max dF/F to find ROIS

b = 5;
siz = size(data_reg);
corr_map = zeros(siz(1),siz(2));
for ix = b:siz(2)-b
    for iy = b:siz(1)-b
        TC = data_reg(iy,ix,:);
        surround = (data_reg(iy-1,ix-1,:)+data_reg(iy-1,ix,:)+data_reg(iy-1,ix+1,:)+data_reg(iy,ix-1,:)+data_reg(iy,ix+1,:)+data_reg(iy+1,ix-1,:)+data_reg(iy+1,ix,:)+data_reg(iy+1,ix+1,:))/8;
        R = corrcoef(TC,surround);
        corr_map(iy,ix) = R(1,2);
    end
end
figure; imagesq(corr_map); colormap(gray)

bwout = imCellEditInteractive(corr_map);
mask_cell = bwlabel(bwout);

%timecourses
data_TC = stackGetTimeCourses(dFoverF_data,mask_cell);
figure; tcOffsetPlot(data_TC)

save('mask&TCDir.mat','mask_cell','data_TC');

%% dF/F (by cell) for each stimulus type

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

% % resp trace for each trial: 8:12+3 frame window
% dFoverF_RspON = zeros(nON+nON+2,nCells,nTrials);
% data_TCplus = data_TC;
% data_TCplus(end+1:end+3,:) = 0;
% for icell = 1:nCells
%     for itrial = 1:nTrials
%         tri = data_TCplus((nON_ind_firsts(itrial)-2):((nON_ind_firsts(itrial)+nON)+2),icell);
%         dFoverF_RspON(:,icell,itrial) = tri;
%     end
% end

%make all positive numbers
% minAdd = abs(min(min(dFoverF_meanRspbyDirbyCell,[],1)));
% dFoverF_meanRspbyDirbyCell = bsxfun(@plus, dFoverF_meanRspbyDirbyCell,minAdd);
% clear minAdd
minAdd = abs(min(min(dFoverF_RspbyCellbyTrial,[],1)));
dFoverF_RspbyCellbyTrial = bsxfun(@plus, dFoverF_RspbyCellbyTrial,minAdd);
clear minAdd
% minAdd = abs(min(min(dFoverF_RspON,[],1)));
% dFoverF_RspON = bsxfun(@plus, dFoverF_RspON,minAdd);
% clear minAdd

%% Preferred orientation

%average response to each orientation
for idir = 1:nStim
    thisDir = Dirs(idir);
    dir_ind = find(DirectionDeg == idir);
    for icell = 1:nCells
        for itrial = 1:nTrials
            

%% Plot overlay of all traces for each direction, for each cell
nCdir = (nDir/2);
    %Dir0 = 0 and 180 deg motion
% for icell = 1:nCells
%     ind = find(DirectionDeg==Dir1(1) | DirectionDeg==Dir2(1));
%     start = 1;
%     for idir = ind
%         Dir0_RspbyCellsbyTrial(:,icell,start) = dFoverF_RspbyCellbyTrial(:,icell,idir);
%         start = start+1;
%     end
% end
%     %Dir45 = 45 and 225 deg motion
% for icell = 1:nCells
%     ind = find(DirectionDeg==Dir1(2) | DirectionDeg==Dir2(2));
%     start = 1;
%     for idir = ind
%         Dir45_RspbyCellsbyTrial(:,icell,start) = dFoverF_RspbyCellbyTrial(:,icell,idir);
%         start = start+1;
%     end
% end
%     %Dir90 = 90 and 270 deg motion
% for icell = 1:nCells
%     ind = find(DirectionDeg==Dir1(3) | DirectionDeg==Dir2(3));
%     start = 1;
%     for idir = ind
%         Dir90_RspbyCellsbyTrial(:,icell,start) = dFoverF_RspbyCellbyTrial(:,icell,idir);
%         start = start+1;
%     end
% end
%     %Dir0 = 135 and 315 deg motion
% for icell = 1:nCells
%     ind = find(DirectionDeg==Dir1(4) | DirectionDeg==Dir2(4));
%     start = 1;
%     for idir = ind
%         Dir135_RspbyCellsbyTrial(:,icell,start) = dFoverF_RspbyCellbyTrial(:,icell,idir);
%         start = start+1;
%     end
% end
% 
% % limit number of cells plotted
% rCells_ind = find(CellDSI > 0.25);
% nCells = size(rCells_ind,2);
% dFoverF_RspbyCellbyTrial_rCells = dFoverF_RspbyCellbyTrial(:,rCells_ind,:);
% 
% for icell = 1:nCells
%     ind = find(DirectionDeg==Dir1(1) | DirectionDeg==Dir2(1));
%     start = 1;
%     for idir = ind
%         Dir0_RspbyCellsbyTrial_rCells(:,icell,start) = dFoverF_RspbyCellbyTrial_rCells(:,icell,idir);
%         start = start+1;
%     end
% end
%     %Dir45 = 45 and 225 deg motion
% for icell = 1:nCells
%     ind = find(DirectionDeg==Dir1(2) | DirectionDeg==Dir2(2));
%     start = 1;
%     for idir = ind
%         Dir45_RspbyCellsbyTrial_rCells(:,icell,start) = dFoverF_RspbyCellbyTrial_rCells(:,icell,idir);
%         start = start+1;
%     end
% end
%     %Dir90 = 90 and 270 deg motion
% for icell = 1:nCells
%     ind = find(DirectionDeg==Dir1(3) | DirectionDeg==Dir2(3));
%     start = 1;
%     for idir = ind
%         Dir90_RspbyCellsbyTrial_rCells(:,icell,start) = dFoverF_RspbyCellbyTrial_rCells(:,icell,idir);
%         start = start+1;
%     end
% end
%     %Dir0 = 135 and 315 deg motion
% for icell = 1:nCells
%     ind = find(DirectionDeg==Dir1(4) | DirectionDeg==Dir2(4));
%     start = 1;
%     for idir = ind
%         Dir135_RspbyCellsbyTrial_rCells(:,icell,start) = dFoverF_RspbyCellbyTrial_rCells(:,icell,idir);
%         start = start+1;
%     end
% end


    %Dir0 = 0 and 180 deg motion
for icell = 1:nCells
    ind = find(DirectionDeg==Dir1(1) | DirectionDeg==Dir2(1));
    start = 1;
    for idir = ind
        Dir0_RspbyCellsbyTrial(:,icell,start) = dFoverF_RspON(:,icell,idir);
        start = start+1;
    end
end
    %Dir45 = 45 and 225 deg motion
for icell = 1:nCells
    ind = find(DirectionDeg==Dir1(2) | DirectionDeg==Dir2(2));
    start = 1;
    for idir = ind
        Dir45_RspbyCellsbyTrial(:,icell,start) = dFoverF_RspON(:,icell,idir);
        start = start+1;
    end
end
    %Dir90 = 90 and 270 deg motion
for icell = 1:nCells
    ind = find(DirectionDeg==Dir1(3) | DirectionDeg==Dir2(3));
    start = 1;
    for idir = ind
        Dir90_RspbyCellsbyTrial(:,icell,start) = dFoverF_RspON(:,icell,idir);
        start = start+1;
    end
end
    %Dir0 = 135 and 315 deg motion
for icell = 1:nCells
    ind = find(DirectionDeg==Dir1(4) | DirectionDeg==Dir2(4));
    start = 1;
    for idir = ind
        Dir135_RspbyCellsbyTrial(:,icell,start) = dFoverF_RspON(:,icell,idir);
        start = start+1;
    end
end

% limit number of cells plotted
rCells_ind = find(CellDSI > 0.25);
nCells = size(rCells_ind,2);
dFoverF_RspON_rCells = dFoverF_RspON(:,rCells_ind,:);

for icell = 1:nCells
    ind = find(DirectionDeg==Dir1(1) | DirectionDeg==Dir2(1));
    start = 1;
    for idir = ind
        Dir0_RspbyCellsbyTrial_rCells(:,icell,start) = dFoverF_RspON_rCells(:,icell,idir);
        start = start+1;
    end
end
    %Dir45 = 45 and 225 deg motion
for icell = 1:nCells
    ind = find(DirectionDeg==Dir1(2) | DirectionDeg==Dir2(2));
    start = 1;
    for idir = ind
        Dir45_RspbyCellsbyTrial_rCells(:,icell,start) = dFoverF_RspON_rCells(:,icell,idir);
        start = start+1;
    end
end
    %Dir90 = 90 and 270 deg motion
for icell = 1:nCells
    ind = find(DirectionDeg==Dir1(3) | DirectionDeg==Dir2(3));
    start = 1;
    for idir = ind
        Dir90_RspbyCellsbyTrial_rCells(:,icell,start) = dFoverF_RspON_rCells(:,icell,idir);
        start = start+1;
    end
end
    %Dir0 = 135 and 315 deg motion
for icell = 1:nCells
    ind = find(DirectionDeg==Dir1(4) | DirectionDeg==Dir2(4));
    start = 1;
    for idir = ind
        Dir135_RspbyCellsbyTrial_rCells(:,icell,start) = dFoverF_RspON_rCells(:,icell,idir);
        start = start+1;
    end
end


% subtract baseline so zero is at nON transition

for icell = 1:nCells
    tri = size(Dir0_RspbyCellsbyTrial_rCells,3);
    for itrial = 1:tri
        Dir0_zeroed(:,icell,itrial) = bsxfun(@minus, Dir0_RspbyCellsbyTrial_rCells(:,icell,itrial), Dir0_RspbyCellsbyTrial_rCells(3,icell,itrial));
    end
end

for icell = 1:nCells
    tri = size(Dir45_RspbyCellsbyTrial_rCells,3);
    for itrial = 1:tri
        Dir45_zeroed(:,icell,itrial) = bsxfun(@minus, Dir45_RspbyCellsbyTrial_rCells(:,icell,itrial), Dir45_RspbyCellsbyTrial_rCells(3,icell,itrial));
    end
end

for icell = 1:nCells
    tri = size(Dir90_RspbyCellsbyTrial_rCells,3);
    for itrial = 1:tri
        Dir90_zeroed(:,icell,itrial) = bsxfun(@minus, Dir90_RspbyCellsbyTrial_rCells(:,icell,itrial), Dir90_RspbyCellsbyTrial_rCells(3,icell,itrial));
    end
end

for icell = 1:nCells
    tri = size(Dir135_RspbyCellsbyTrial_rCells,3);
    for itrial = 1:tri
        Dir135_zeroed(:,icell,itrial) = bsxfun(@minus, Dir135_RspbyCellsbyTrial_rCells(:,icell,itrial), Dir135_RspbyCellsbyTrial_rCells(3,icell,itrial));
    end
end

%plot
figure;
for icell = 1:21
    subplot(5,5,icell)
    x = squeeze(Dir0_zeroed(:,icell,:));
    plot(x)
end
figure;
for icell = 1:21
    subplot(5,5,icell)
    x = squeeze(Dir45_zeroed(:,icell,:));
    plot(x)
end             
figure;
for icell = 1:21
    subplot(5,5,icell)
    x = squeeze(Dir90_zeroed(:,icell,:));
    plot(x)
end
figure;
for icell = 1:21
    subplot(5,5,icell)
    x = squeeze(Dir135_zeroed(:,icell,:));
    plot(x)
end

%% Direction preferences

%prefered orientation
CellPref_val = zeros(1,nCells);
CellPref_ind = zeros(1,nCells);
for icell = 1:nCells
    [CellPref_val(:,icell) CellPref_ind(:,icell)] = max(dFoverF_meanRspbyDirbyCell(:,icell),[],1);
end

CellPref_deg = Dir(CellPref_ind);

%OSI
CellOrth_deg = zeros(size(CellPref_deg));
for icell = 1:nCells
    if CellPref_deg(icell) < 180
        CellOrth_deg(icell) = CellPref_deg(icell) + 180;
    elseif CellPref_deg(icell) >= 180
        CellOrth_deg(icell) = CellPref_deg(icell)-180;
    end
end
for icell = 1:nCells
    CellOrth_ind(icell) = find(Dir == CellOrth_deg(icell));
end

CellRspPref = zeros(1,nCells);
for icell = 1:nCells
    dir_ind = CellPref_ind(icell);
    CellRspPref(:,icell) = dFoverF_meanRspbyDirbyCell(dir_ind,icell);
end

CellRspOrth = zeros(1,nCells);
for icell = 1:nCells
    dir_ind = CellOrth_ind(icell);
    CellRspOrth(:,icell) = dFoverF_meanRspbyDirbyCell(dir_ind,icell);
end

CellOSI = (CellRspPref - CellRspOrth)./(CellRspPref + CellRspOrth);

%preferred direction is same as preferred orientation

%DSI
CellMinus90plus_deg = zeros(size(CellPref_deg));
CellMinus90minus_deg = zeros(size(CellPref_deg));
for icell = 1:nCells
    if CellPref_deg(icell) >= 90 & CellPref_deg(icell) < 270
        CellMinus90plus_deg(icell) = CellPref_deg(icell) + 90;
        CellMinus90minus_deg(icell) = CellPref_deg(icell)-90;
    elseif CellPref_deg(icell) < 90
        CellMinus90plus_deg(icell) = CellPref_deg(icell)+ 90;
        CellMinus90minus_deg(icell) = (CellPref_deg(icell) - 90) +360;
    elseif CellPref_deg(icell) >= 270
        CellMinus90plus_deg(icell) = (CellPref_deg(icell)+ 90) -360;
        CellMinus90minus_deg(icell) = CellPref_deg(icell) - 90;
    end
end

for icell = 1:nCells
    CellMinus90plus_ind(icell) = find(Dir == CellMinus90plus_deg(icell));
    CellMinus90minus_ind(icell) = find(Dir == CellMinus90minus_deg(icell));
end


CellRspMinus90plus = zeros(1,nCells);
CellRspMinus90minus = zeros(1,nCells);
for icell = 1:nCells
    dirPlus_ind = CellMinus90plus_ind(icell);
    dirMinus_ind = CellMinus90minus_ind(icell);
    CellRspMinus90plus(:,icell) = dFoverF_meanRspbyDirbyCell(dirPlus_ind,icell);
    CellRspMinus90minus(:,icell) = dFoverF_meanRspbyDirbyCell(dirMinus_ind,icell);
    CellRspMinus90(:,icell) = (CellRspMinus90plus(:,icell) + CellRspMinus90minus(:,icell))/2;
end

CellDSI = (CellRspPref - CellRspMinus90)./(CellRspPref + CellRspMinus90);

%view tuning for all cells
figure;
start = 1;
for icell = 51:75
    subplot(5,5,start);
    plot(dFoverF_meanRspbyDirbyCell(:,icell), 'r');
    hold on;
    axis([0 13 -0.1 0.5]);
    hline(0,'k-');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end

figure;
start = 1;
for icell = 26:nCells
    subplot(5,5,start);
    plot(dFoverF_meanRspbyDirbyCell(:,icell), 'r');
    hold on;
    axis([0 13 -0.1 0.5]);
    hline(0,'k-');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end

%cells tuned to flashing stimulus
deg0_ind = find(CellPref_deg == 0);
Rsp_deg0 = dFoverF_meanRspbyDirbyCell(:,deg0_ind);
figure;
start = 1;
siz = size(Rsp_deg0,2);
for icell = deg0_ind
    for ipref = 1:siz
    subplot(3,3,start);
    plot(Rsp_deg0(:,ipref), 'r');
    hold on;
    axis([0 13 -0.1 0.5]);
    hline(0,'k-');
    title(['Cell ' num2str(icell)]);
    start = start+1;
    end
end

%cells tuned to range of change stimulus
degLessThan0_ind = find(CellPref_deg>250 & CellPref_deg<360);
Rsp_degLessThan0 = dFoverF_meanRspbyDirbyCell(:,degLessThan0_ind);
figure;
start = 1;
siz = size(Rsp_degLessThan0,2);
for icell = degLessThan0_ind
        subplot(3,3,start);
        plot(Rsp_degLessThan0(:,start), 'r');
        hold on;
        axis([0 13 -0.1 0.5]);
        hline(0,'k-');
        title(['Cell ' num2str(icell)]);
        start = start+1;
end


