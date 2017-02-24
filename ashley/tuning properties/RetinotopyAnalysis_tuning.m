rc = behavConstsAV;
awTempRet
iexp = 2;

% for iexp = 1:size(expt,2);
SubNum = expt(iexp).SubNum;
expDate = expt(iexp).date;
expTime = expt(iexp).rettuning{2,:};
retFolder = expt(iexp).rettuning{1,:};
mouse = expt(iexp).mouse;
fName = [retFolder '_000_000'];

[input] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,retFolder,fName);    

fnout = fullfile(rc.ashleyAnalysis,mouse,expt(iexp).folder,expDate,retFolder);

%%
down = 10;

if exist(fullfile(fnout,'reg_img'))
    data_reg = readtiff(fullfile(fnout,'retTuning.tif'));
else
    [input, data] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,retFolder,fName);

    data_down = stackGroupProject(data,down);
    clear data

    data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
    data = data_sub;
    clear data_sub

    data_avg = mean(data(:,:,200:210),3);
    figure; imagesq(data_avg); colormap(gray)

    [out data_reg] = stackRegister(data, data_avg);
    clear data
    writetiff(data_reg,fullfile(fnout,'retTuning.tif'));
    save(fullfile(fnout,'reg_img'),'data_avg');
end

try
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', expDate, retFolder);
    cd(filedir);
catch
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging');
    cd(filedir)
    mkdir(expDate,retFolder)
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', expDate, retFolder);
    cd(filedir);
end

%% 
fr_1s = expt(iexp).frame_rate./down;
on = double(input.nScansOn/down);
off = double(input.nScansOff/down);
base_win = off-fr_1s+1:off;
resp_win = off+1:off+(2.*fr_1s);

nfr_tr = on+off;
nfr_session = floor(size(data_reg,3)./nfr_tr).*nfr_tr;
data_reg = data_reg(:,:,1:nfr_session);

trStart = 1:on+off:nfr_session;

off_ind = linspaceNDim(trStart,trStart+off-1,off);

nTrials = input.trialSinceReset;
nAz = double(input.gratingAzimuthStepN);
nEl = double(input.gratingElevationStepN);
tAzs = cell2mat_padded(input.tGratingAzimuthDeg);
tEls = cell2mat_padded(input.tGratingElevationDeg);
[azind azs] = findgroups(tAzs);
[elind els] = findgroups(tEls);

nPos = nAz*nEl;
tPos = zeros(1,nTrials);
az_mat = zeros(nPos,1);
el_mat = zeros(nPos,1);
for iaz = 1:nAz
    a_ind = find(azind == iaz);
    for iel = 1:nEl
        e_ind = find(elind == iel);
        pos_ind = iaz + (nAz * (iel-1));
        tPos(intersect(a_ind,e_ind)) = pos_ind;
        az_mat(pos_ind) = azs(iaz);
        el_mat(pos_ind) = els(iel);
    end
end

% % % % %%
% % % % nRep = size(data_reg,3)./((nON+nOFF)*nStim);
% % % % nTrials = (nStim.*nRep);
% % % % if (mod(nRep,1)) >0
% % % %     nframes = floor(nRep)*((nON+nOFF)*nStim)
% % % %     data_reg = data_reg(:,:,1:nframes);
% % % %     nRep = size(data_reg,3)./((nON+nOFF)*nStim);
% % % %     nTrials = (nStim.*nRep);
% % % % end
% % % % 

if exist(fullfile(fnout,'ret_mask.mat'))
    load(fullfile(fnout,'ret_mask.mat'))
else
    Fimg = mean(data_reg(:,:,off_ind(:)),3);
    dFimg = bsxfun(@minus,data_reg,Fimg);
    dFFimg = bsxfun(@rdivide,data_reg,Fimg);
    maxDFFimg = max(dFFimg,[],3);

    bwout = imCellEditInteractive(maxDFFimg);
    mask_cell = bwlabel(bwout);
    save(fullfile(fnout,'ret_mask'),'mask_cell')
end
%%

% % % %write tifs for sorted frames
% % % VSR = 1;
% % % run('sortTrialsAvg_writetiffs.m')

%% get timecourses and subtract neuropil
% get timecourses
% % % try
% % % dataTC = stackGetTimeCourses(data_reg, mask_cell);
% % % catch
% % %     dataTC = stackGetTimeCourses(data_reg, mask_boutons);
% % % end

dataTC = stackGetTimeCourses(data_reg, mask_cell);

% get neuropil timecourses
buf = 4;
np = 6;
nCells = size(dataTC,2);
neuropil = imCellNeuropil(mask_cell,buf,np);


npTC = zeros(size(dataTC));
for i = 1:nCells
    tempNPmask = squeeze(neuropil(:,:,i));
    if sum(sum(tempNPmask)) > 0
    npTC(:,i) = stackGetTimeCourses(data_reg,tempNPmask);
    end
end

% npTC = stackGetTimeCourses(data_reg,neuropil);


dataTC_mavg = tsmovavg(dataTC,'s',10,1);
npTC_mavg = tsmovavg(npTC,'s',10,1);

% down = 10;
% dataTC_down = reshape(mean(reshape(dataTC',size(dataTC',1),down,size(dataTC',2)/down),2),size(dataTC',1),size(dataTC',2)/down)';
% npTC_down = reshape(mean(reshape(npTC',size(npTC',1),down,size(npTC',2)/down),2),size(npTC',1),size(npTC',2)/down)';

ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(dataTC_mavg-tcRemoveDC(npTC_mavg*ii(i)));
end
[max_skew ind] =  max(x,[],1);
% skew(buf,:) = max_skew;
np_w = 0.01*ind;
dataTCsub = dataTC-bsxfun(@times,tcRemoveDC(npTC),np_w);
data_TC = dataTCsub;
clear data_reg

% % % %% vis stim parameters
% % % VSsize = input.gratingDiameterDeg;
% % % VSdirectionDeg = input.gratingDirectionDeg;
% % % tAz = double(cell2mat(input.tGratingAzimuthDeg));
% % % Az = unique(tAz);
% % % [nAz azPos] = histc(tAz,Az);
% % % tEl = double(cell2mat(input.tGratingElevationDeg));
% % % El = unique(tEl);
% % % [nEl elPos] = histc(tEl,El);
% % % if any(tEl == 0) | any(tAz == 0)
% % %     if any(tEl == 0)
% % %         tEl2 = tEl+1;
% % %     else
% % %         tEl2 = tEl;
% % %     end
% % %     if any(tAz == 0)
% % %         tAz2 = tAz+1;
% % %     else
% % %         tAz2 = tAz;
% % %     end
% % %     pos = cart2pol(tAz2,tEl2);
% % % else
% % %     pos = cart2pol(tAz,tEl);
% % % end
% % % [n posN] = histc(pos,unique(pos));
% % % Rets = unique(pos);
% % % 
% % % AzPos = NaN(1,nStim);
% % % ElPos = NaN(1,nStim);
% % % for i = 1:nStim
% % %     AzPos(i) = unique(tAz(pos == Rets(i)));
% % %     ElPos(i) = unique(tEl(pos == Rets(i)));
% % % end
% % % 


%% dF/F by trial
data_tr = reshape(data_TC,nfr_tr,nTrials,nCells);

F = mean(data_tr(base_win,:,:),1);
dF = bsxfun(@minus,data_tr,F);
dFF = bsxfun(@rdivide,dF,F);

%% stim responses
resp = zeros(nfr_tr,nCells,nPos);
resperr = zeros(nfr_tr,nCells,nPos);
for ipos = 1:nPos
    ind = tPos == ipos;
    resp(:,:,ipos) = squeeze(mean(dFF(:,ind,:),2));    
end

posR = squeeze(mean(resp(resp_win,:,:),1));
posRerr = squeeze(mean(resperr(resp_win,:,:),1));

meanR = mean(posR,1);
Rimg = reshape(meanR',nAz,nEl)';

trS = chop((-off+1:on)./double(expt(iexp).frame_rate/down),2);
setFigParams4Print('landscape')
figure;
for ipos = 1:nPos
   subplot(nEl,nAz,double(ipos))
   plot(trS,resp(:,:,ipos))
   xlim([trS(1) trS(end)])
   ylim([-0.5 2])
   title([num2str(az_mat(ipos)) '/' num2str(el_mat(ipos)) ' az/el'])
   xlabel('s')
   ylabel('dF/F')
end
print(fullfile(fnout,'cell_ret_resp'),'-dpdf','-fillpage');

setFigParams4Print('portrait')
figure;
h = imagesc(Rimg);
h.Parent.XTick = 1:nAz;
h.Parent.YTick = 1:nEl;
h.Parent.XTickLabel = azs;
h.Parent.YTickLabel = els;
colorbar
caxis([0 max(Rimg(:))+0.05])
title('mean resp all cells')
print(fullfile(fnout,'map_ret_resp'),'-dpdf','-fillpage');


% % % 
% % % %% sort data by trial type
% % % dFoverFCellsTrials = zeros(10+on,size(dFoverF_data,2),nTrials);
% % % for i = 1:nTrials
% % %     dFoverFCellsTrials(:,:,i) = dFoverF_data(stimON_ind(i)-10:stimON_ind(i)+(on-1),:);
% % % end
% % % 
% % % dFoverF_meanRetResp = zeros(size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2),nStim);
% % % errbarRets = zeros(size(data_TC,2),nStim);
% % % for i = 1:nStim
% % %     trials = find(pos(:,1:nTrials) == Rets(i));
% % %     dFoverF_meanRetResp(:,:,i) = mean(dFoverFCellsTrials(:,:,trials),3);
% % %     errbarRets(:,i) = std(mean(dFoverFCellsTrials(11:16,:,trials),1),[],3)/sqrt(size(dFoverFCellsTrials(11:16,:,trials),3));
% % % end
% % % 
% % % figure;
% % % for i = 1:nStim
% % %     plot(dFoverF_meanRetResp(:,10,i));
% % %     hold on
% % % end
% % % 
% % % %% plot tuning curves
% % % dFoverF_meanOFFRetResp = (squeeze(mean(dFoverF_meanRetResp(1:10,:,:),1)));
% % % 
% % % RetRespPerCell = (squeeze(mean(dFoverF_meanRetResp(11:end,:,:),1)));
% % % 
% % % % az x el matrix per cell
% % % resp = squeeze(mean(dFoverFCellsTrials(11:end,:,:),1));
% % % resp_AzElMat = zeros(size(resp,1),nRep,length(El),length(Az));
% % % 
% % % for iaz = 1:length(Az)
% % %     for iel = 1:length(El)
% % %         resp_AzElMat(:,:,iel,iaz) = resp(:,find(tEl == El(iel) & tAz == Az(iaz)));
% % %     end
% % % end
% % % 
% % % meanResp_AzElMat = squeeze(mean(resp_AzElMat,2));
% % % errResp_AzElMat = squeeze(std(resp_AzElMat,[],2)/sqrt(double(nRep)));
% % % 
% % % %%
% % % figure;
% % % imagesc(squeeze(mean(meanResp_AzElMat,1)));
% % % set(gca,'xTick', [1:length(Az)])
% % % set(gca,'xTickLabel', cellfun(@num2str,num2cell(Az),'UniformOutput',false))
% % % set(gca,'yTick', [1:length(El)])
% % % set(gca,'yTickLabel', cellfun(@num2str,num2cell(El),'UniformOutput',false))
% % % xlabel('Az')
% % % ylabel('El')
% % % colorbar
% % % 
% % % runstr = expt(iexp).runs(1,:);
% % % if expt(iexp).nrun>1
% % %     for irun = 2:expt(iexp).nrun
% % %         runstr = [runstr '-' expt(iexp).runs(irun,:)];
% % %     end
% % % end
% % % var = whos('-file',fullfile('Z:\analysis\',mouse,'two-photon imaging', expDate, runstr,[mouse '-' expDate '-' runstr '-comboInputDataTCplusVar.mat']));
% % % load(fullfile('Z:\analysis\',mouse,'two-photon imaging', expDate, runstr,[mouse '-' expDate '-' runstr '-comboInputDataTCplusVar.mat']), var(structfind(var,'name','input')).name)
% % % retUsedAz = input.gratingAzimuthDeg;
% % % retUsedEl = input.gratingElevationDeg;
% % % 
% % % title({'Retinotopy';[SubNum '-' expDate]; ['(Az,El) = (' num2str(retUsedAz) ',' num2str(retUsedEl) ')'] })
% % % 
% % % set(0,'defaultfigurepaperorientation','portrait');
% % % set(0,'defaultfigurepapersize',[8.5 11]);
% % % set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
% % % 
% % % % get position used for behavior
% % % 
% % % 
% % % print(fullfile('Z:\analysis\',mouse,'two-photon imaging', expDate, retFolder,'RetPreferences'), '-dpdf')
%% save tuning info
save(fullfile(fnout,'Timecourses.mat') ,'data_TC','npTC')
save(fullfile(fnout,'RetResp.mat'),'resp','azs','els','tPos')

%%