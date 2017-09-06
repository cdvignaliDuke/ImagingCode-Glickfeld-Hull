clear all
close all
ds = '_V1gad';
rc = behavConstsAV;
eval(['awFSAVdatasets' ds])
slct_expt = 7;
%%
for iexp = slct_expt
    if expt(iexp).greenredsimultaneous ~= 1
        continue
    end

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
dirFolder = expt(iexp).dirtuning;
dirTime = expt(iexp).dirtuning_time;
fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
%% load all behavior data
% concatenate data files
for irun = 1:expt(iexp).nrun
 
    runFolder = expt(iexp).runs(irun,:);
    expTime = expt(iexp).time_mat(irun,:);
    fName = [runFolder '_000_000'];
    if expt(iexp).redChannelOn
%         greendata_temp = loadsbx_choosepmt(1,mouse,expDate,runFolder,fName);
        data_temp = loadsbx_choosepmt(2,mouse,expDate,runFolder,fName);
        input = loadMworksFile(SubNum,expDate,expTime);
    else
        continue
    end
    
    if irun == 1
        data_bx = data_temp;
        input_bx = input;
    else
        data_bx = cat(3,data_bx,data_temp);
        input_bx = [input_bx input];
    end
    clear data_temp input
end
input_bx = concatenateDataBlocks(input_bx);
nfr_bx = size(data_bx,3);

% remove negative data by subtraction
data_sub = data_bx-min(min(min(data_bx,[],1),[],2),[],3);
data_bx = data_sub;
clear data_sub

% load tuning data
down = 10;
fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,dirFolder);
fName = [dirFolder '_000_000'];

data_tun = loadsbx_choosepmt(2,mouse,expDate,dirFolder,fName);
input_tun = loadMworksFile(SubNum,expDate,dirTime);
% down-sample
data_tun_down = stackGroupProject(data_tun,down);
clear data_tun

% remove negative data by subtraction
data_sub = data_tun_down-min(min(min(data_tun_down,[],1),[],2),[],3);
data_tun_down = data_sub;
clear data_sub

nfr_tun = size(data_tun_down,3);

figure
randframes = 1:100;
for i = 1:9
    randframes = randframes+100;
    if randframes(end) > nfr_bx
        continue
    end
    img_temp = mean(data_bx(:,:,randframes),3);
    subplot(3,3,i)
    imagesc(img_temp)
    title(sprintf('frames %s:%s',num2str(randframes(1)),num2str(randframes(end))))
end

frames4Reg = 201:300;%input('Enter which frames to use for registration in the form x:x+100: ');

image4Reg = mean(data_bx(:,:,frames4Reg),3);

[out_bx, data_bx_reg] = stackRegister(data_bx,image4Reg);
clear data_bx

[out_tun, data_tun_reg] = stackRegister(data_tun_down,image4Reg);
clear data_tun
% 
% fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
if ~exist(fnout,'dir')
    mkdir(fnout)
end
save(fullfile(fnout,'reg2RedOuts&Img.mat'),'out_bx','out_tun', 'image4Reg','nfr_bx','nfr_tun')
%
end