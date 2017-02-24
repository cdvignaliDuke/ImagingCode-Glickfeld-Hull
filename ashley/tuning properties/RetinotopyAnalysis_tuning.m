rc = behavConstsAV;
awTempRet
iexp = 1;

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

%% get timecourses and subtract neuropil

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

%% save tuning info
save(fullfile(fnout,'Timecourses.mat') ,'data_TC','npTC')
save(fullfile(fnout,'RetResp.mat'),'resp','azs','els','tPos')

%%