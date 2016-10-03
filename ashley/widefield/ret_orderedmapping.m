mouse = 'AW58';
SubNum = '658';
expdate = '160831';
foldername = 'msFRONT_ret6pos_1';
filename = 'msFRONT_ret6pos_1.tif';

nON = 10;
nOFF = nON;
nStim = 6;
nframes = 1200;
ntrials = nframes/(nON+nOFF);
nRep = ntrials/nStim;
trans_win = nOFF+1:nOFF+5;

fnin = fullfile('Z:\data',mouse,'widefield imaging',[mouse '_' expdate]);
fnout = fullfile('Z:\Analysis',mouse,'widefield imaging',expdate);
mkdir(fnout);

data = readtiff(fullfile(fnin,foldername,filename));
data = double(data);

dFoverF = zeros(size(data,1),size(data,2),nOFF+nOFF,ntrials);
for i = 1:ntrials
    prevFrames = (nON+nOFF)*(i-1)
    F = mean(data(:,:,prevFrames+nOFF-4:prevFrames+nOFF),3);
    dFoverF(:,:,:,i) = bsxfun(@rdivide, bsxfun(@minus,data(:,:,prevFrames+1:prevFrames+nON+nOFF),F), F);
end

trial_ind = zeros(nStim,nRep);
for i = 1:nStim
    trial_ind(i,:) = i:nStim:ntrials;
end   

dFoverF_avg_trial = zeros(size(data,1),size(data,2),nStim);
for i = 1:nStim
    dFoverF_avg_trial(:,:,i) = mean(mean(dFoverF(:,:,trans_win,trial_ind(i,:)),4),3);
end

figure;
for i = 1:nStim
    subplot(2,3,i)
    imagesc(dFoverF_avg_trial(:,:,i))
end

