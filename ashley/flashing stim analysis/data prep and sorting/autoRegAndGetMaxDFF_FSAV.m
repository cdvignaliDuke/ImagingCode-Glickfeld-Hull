clear all
close all
ds = '_naive100ms';
rc = behavConstsAV;
eval(['awFSAVdatasets' ds])
slct_expt = 1:size(expt,2);
%%
for iexp = slct_expt

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
dirFolder = expt(iexp).dirtuning;
dirTime = expt(iexp).dirtuning_time;

%% load all behavior data
% concatenate data files
for irun = 1:expt(iexp).nrun

    if strcmp(ds, '_V1') & irun == 3 & strcmp(expDate,'150508')
        continue
    elseif strcmp(ds, '_V1') & irun == 2 & strcmp(expDate, '150528')
        continue
    end
    
    runFolder = expt(iexp).runs(irun,:);
    expTime = expt(iexp).time_mat(irun,:);
    fName = [runFolder '_000_000'];
    [input, data_temp, t] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,runFolder,fName);
    disp(t)
    
    if irun == 1
        data_bx = data_temp;
        input_bx = input;
    else
        data_bx = cat(3,data_bx,data_temp);
        try
            input_bx = [input_bx input];
        catch
            inpNames1 = fieldnames(input_bx);
            inpNames2 = fieldnames(input);
            inpLong = gt(length(inpNames1),length(inpNames2));
            if inpLong == 1
                inpPlusInd = ismember(inpNames1,inpNames2);
                inpPlus = inpNames1(~inpPlusInd);
                for i = 1:length(inpPlus)
                    input.(genvarname(inpPlus{i})) = cell(1,input.trialSinceReset);
                end
            else
                inpPlusInd = ismember(inpNames2,inpNames1);
                inpPlus = inpNames2(~inpPlusInd);
                for i = 1:length(inpPlus)
                    input_bx.(char(genvarname(inpPlus(i)))) = cell(1,length(inpNames1));
                end
            end
            input_temp = [input_bx input];
        end
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

if any(strcmp(fieldnames(expt),'nTunFrames'))
    [input_tun, data_tun] = Load_SBXdataPlusMWorksData(SubNum,expDate,dirTime,mouse,dirFolder,fName,expt(iexp).nTunFrames);  
else
    [input_tun, data_tun] = Load_SBXdataPlusMWorksData(SubNum,expDate,dirTime,mouse,dirFolder,fName);
end

% down-sample
data_tun_down = stackGroupProject(data_tun,down);
clear data_tun

% remove negative data by subtraction
data_sub = data_tun_down-min(min(min(data_tun_down,[],1),[],2),[],3);
data_tun_down = data_sub;
clear data_sub

nfr_tun = size(data_tun_down,3);

% get average frames
d = 100;
nd = floor(size(data_bx,3)/d);
xpix = size(data_bx,2);
ypix = size(data_bx,1);

data_d = zeros(ypix,xpix,nd);
for id = 1:nd
   data_d(:,:,id) = mean(data_bx(:,:,(d*(id-1))+1:id*d),3); 
end

data_d = cat(3, data_d, data_tun_down);

% find most correlated frame in unregistered data
data_d_2D = reshape(data_d,xpix*ypix,nd+nfr_tun);
data_corr = corrcoef(data_d_2D);
clear data_d_2D 

[max_val max_corr_ind] = max(mean(data_corr));
clear data_corr
% register downsampled data
[out data_d_reg] = stackRegister(data_d,data_d(:,:,max_corr_ind));

data_reg_2D = reshape(data_d_reg,xpix*ypix,nd+nfr_tun);
data_corr2 = corrcoef(data_reg_2D);
% find max correlation frame in registered data
[max_corr2_val max_corr2_ind] = max(mean(data_corr2));

data_corr_img = data_d_reg(:,:,max_corr2_ind);

clear data_d data_d_reg data_reg_2D data_corr2 
% register all behavior data and save outs and reg image
[out_bx data_bx_reg] = stackRegister(data_bx,data_corr_img);
clear data_bx

[out_tun data_tun_reg] = stackRegister(data_tun_down,data_corr_img);
clear data_tun

fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
if ~exist(fnout,'dir')
    mkdir(fnout)
end
save(fullfile(fnout,'regOuts&Img.mat'),'out_bx','out_tun', 'data_corr_img','nfr_bx','nfr_tun')
%% find motion trials
% tuning
motion_fig = figure; setFigParams4Print('portrait')
subplot(2,2,1)
histogram(out_tun(:,1))
xlabel('corr coef')
subplot(2,2,2)
histogram(out_tun(:,2))
xlabel('phase diff')
subplot(2,2,3)
histogram(out_tun(:,3))
xlabel('row shift')
subplot(2,2,4)
histogram(out_tun(:,4))
xlabel('col shift')

tMotion_tun = cat(1,find(abs(out_tun(:,3)) > 10), find(abs(out_tun(:,4)) > 10));

nMotion = length(tMotion_tun);
pMotion = num2str(round((nMotion/nfr_tun)*100));

suptitle({[mouse '-' expDate '-dir tuning'];[num2str(nMotion) ' (' pMotion '%) frames with large translation']})

if ~exist(fullfile(rc.ashleyAnalysis,'FSAV Summaries',['awFSAVdatasets' ds],'motion histograms'),'dir')
    mkdir(fullfile(rc.ashleyAnalysis,'FSAV Summaries',['awFSAVdatasets' ds],'motion histograms'))
end
print(fullfile(rc.ashleyAnalysis,'FSAV Summaries',['awFSAVdatasets' ds],'motion histograms',[mouse '_' expDate '_tun']),'-dpdf')

% behavior
motion_fig = figure; setFigParams4Print('portrait')
subplot(2,2,1)
histogram(out_bx(:,1))
xlabel('corr coef')
subplot(2,2,2)
histogram(out_bx(:,2))
xlabel('phase diff')
subplot(2,2,3)
histogram(out_bx(:,3))
xlabel('row shift')
subplot(2,2,4)
histogram(out_bx(:,4))
xlabel('col shift')

tMotion_bx = cat(1,find(abs(out_bx(:,3)) > 10), find(abs(out_bx(:,4)) > 10));

nMotion = length(tMotion_bx);
pMotion = num2str(round((nMotion/nfr_bx)*100));

suptitle({[mouse '-' expDate '-behavior'];[num2str(nMotion) ' (' pMotion '%) frames with large translation']})

fnout = fullfile(rc.ashleyAnalysis,'FSAV Summaries',['awFSAVdatasets' ds]);

print(fullfile(fnout,'motion histograms',[mouse '_' expDate '_bx']),'-dpdf')
%% break up data, set nans for motion frames;

data_bx_reg(:,:,tMotion_bx) = NaN;
data_bx_tun(:,:,tMotion_tun) = NaN;


%% ********* behavior max dF/F ************

fnbx = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
% get behavior max dF/F
taskMaxDFF

%% ********* tuning max dF/F ************
%%
fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,dirFolder);

dirTuningMaxDFF
end