clear all
close all
ds = 'awData_audMod_V13trialtypes';
slct_exp = 2;
regDownSampN = 100; % should be calibrated depending on downsample rate
%%
rc = behavConstsAV;
eval(ds)
for iexp = slct_exp

subnum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
down = expt(iexp).downSampleRate;
for irun = 1:expt(iexp).nrun
    
    runFolder = expt(iexp).runs(irun,:);
    expTime = expt(iexp).time_mat(irun,:);
    fName = [runFolder '_000_000'];
    fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,runFolder);
    if ~exist(fnout,'dir')
        mkdir(fnout)
    end
    [input, data] = Load_SBXdataPlusMWorksData(subnum,expDate,expTime,mouse,runFolder,fName);
    
    % down-sample
    if down > 1
        data_down = stackGroupProject(data,down);
        clear data
    else
        data_down = data;
        clear data
    end
    % remove negative data by subtraction
    data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
    clear data_down

    % get average frames
    d = regDownSampN;
    nd = floor(size(data_sub,3)/d);
    xpix = size(data_sub,2);
    ypix = size(data_sub,1);
    nfr = size(data_sub,3);


    data_d = zeros(ypix,xpix,nd);
    for id = 1:nd
       data_d(:,:,id) = mean(data_sub(:,:,(d*(id-1))+1:id*d),3); 
    end
    
    if irun > 1
        load(fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate,expt(iexp).regImg,'regOuts&Img.mat'))
        
        [out_reg data_reg] = stackRegister(data_sub,data_corr_img);
        clear data_sub

        save(fullfile(fnout,'regOuts&Img.mat'),'out_reg','data_corr_img')        
    else
    % find most correlated frame in unregistered data
    data_d_2D = reshape(data_d,xpix*ypix,nd);
    data_corr = corrcoef(data_d_2D);
    clear data_d_2D 

    [max_val max_corr_ind] = max(mean(data_corr));
    clear data_corr
    
    % register downsampled data
    [out data_d_reg] = stackRegister(data_d,data_d(:,:,max_corr_ind));

    % find max correlation frame in registered data
    data_reg_2D = reshape(data_d_reg,xpix*ypix,nd);
    data_corr2 = corrcoef(data_reg_2D);
    [max_corr2_val max_corr2_ind] = max(mean(data_corr2));

    data_corr_img = data_d_reg(:,:,max_corr2_ind);
    clear data_d data_d_reg data_reg_2D data_corr2
    
    % register raw data and save outs and reg image
    [out_reg data_reg] = stackRegister(data_sub,data_corr_img);
    clear data_sub
    
    save(fullfile(fnout,'regOuts&Img.mat'),'out_reg','data_corr_img')
    end

% % %% find motion trials
% % % tuning
% % motion_fig = figure; setFigParams4Print('portrait')
% % subplot(2,2,1)
% % histogram(out_tun(:,1))
% % xlabel('corr coef')
% % subplot(2,2,2)
% % histogram(out_tun(:,2))
% % xlabel('phase diff')
% % subplot(2,2,3)
% % histogram(out_tun(:,3))
% % xlabel('row shift')
% % subplot(2,2,4)
% % histogram(out_tun(:,4))
% % xlabel('col shift')
% % 
% % tMotion_tun = cat(1,find(abs(out_tun(:,3)) > 10), find(abs(out_tun(:,4)) > 10));
% % 
% % nMotion = length(tMotion_tun);
% % pMotion = num2str(round((nMotion/nfr_tun)*100));
% % 
% % suptitle({[mouse '-' expDate '-dir tuning'];[num2str(nMotion) ' (' pMotion '%) frames with large translation']})
% % 
% % if ~exist(fullfile(rc.ashleyAnalysis,'FSAV Summaries',['awFSAVdatasets' ds],'motion histograms'),'dir')
% %     mkdir(fullfile(rc.ashleyAnalysis,'FSAV Summaries',['awFSAVdatasets' ds],'motion histograms'))
% % end
% % print(fullfile(rc.ashleyAnalysis,'FSAV Summaries',['awFSAVdatasets' ds],'motion histograms',[mouse '_' expDate '_tun']),'-dpdf')
% % %% break up data, set nans for motion frames;
% % 
% % data_reg(:,:,tMotion_bx) = NaN;


%% ********* tuning max dF/F ************
data_reg = double(data_reg);
if expt(iexp).tunVar == 0 & expt(iexp).exptAlignVar == 0
    maxDFF = quickMaxDFF(data_reg);
    %% figure
    figure; setFigParams4Print('portrait')
    colormap gray
    imagesc(maxDFF)

    if ~exist(fullfile(fnout, 'max images'),'dir')
        mkdir(fullfile(fnout, 'max images'))
    end
    print(fullfile(fnout, 'max images'),'-dpdf')
    writetiff(maxDFF, fullfile(fnout,'max images'))
elseif expt(iexp).tunVar ~= 0
    tunVar = expt(iexp).tunVar;
    maxDFF = quickTunMaxDFF(data_reg,input,down,tunVar,fnout,1);
elseif expt(iexp).exptAlignVar ~= 0
    exptVar = expt(iexp).exptAlignVar;
    maxDFF = quickExptMaxDFF(data_reg,input,down,exptVar,expt(iexp).frame_rate,fnout,1);
end
clear data_reg
end
end
