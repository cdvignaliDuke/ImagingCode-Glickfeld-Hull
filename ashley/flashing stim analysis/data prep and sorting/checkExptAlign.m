clear all
close all
fnout = 'Z:\Analysis\FSAV Summaries\awFSAVdatasets_V1\check alignment';
rc = behavConstsAV;
awFSAVdatasets_V1

for iexp = 2:4
    
    SubNum = expt(iexp).SubNum;
    mouse = expt(iexp).mouse;
    expDate = expt(iexp).date;
    dirFolder = expt(iexp).dirtuning;
    dirTime = expt(iexp).dirtuning_time;
    down = 10;
    
    %% load alignment
    fnalign = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);

    load(fullfile(fnalign,'regOuts&Img.mat'))
    
    %% load 100 frames from each behavior dataset and tuning dataset
    
    for irun = 1:expt(iexp).nrun

        if irun == 3 & strcmp(expDate,'150508')
            continue
        end

        runFolder = expt(iexp).runs(irun,:);
        expTime = expt(iexp).time_mat(irun,:);
        fName = [runFolder '_000_000'];
        [input, data_temp] = Load_SBXdataPlusMWorksData(SubNum,expDate,expTime,mouse,runFolder,fName,100);  
        
        load([fName '.mat'])
        nfr = info.config.frames;
        if irun == 1
            fr_ind = 1:100;
        else
            fr_ind = fr_ind+nfr;
        end
        
        out_temp = out_bx(fr_ind,:);
        [out data_temp] = stackRegister_MA(data_temp,[],[],out_temp);
        
        data = stackGroupProject(data_temp,down);
        clear data_temp
        
        data = mean(data,3);
        if irun == 1
            data_means = zeros(size(data,1),size(data,2),expt(iexp).nrun+1);
            data_mean(:,:,irun) = data;
        else
            data_mean(:,:,irun) = data;
        end
        clear input
    end
    
    fName = [dirFolder '_000_000'];
    [input_tun, data_tun] = Load_SBXdataPlusMWorksData(SubNum,expDate,dirTime,mouse,dirFolder,fName,100);  
    
    data = stackGroupProject(data_tun,down);
    clear data_tun input_tun
    
    out_temp = out_tun(1:10,:);
    [out data] = stackRegister_MA(data,[],[],out_temp);
    
    data = mean(data,3);
    
    data_mean(:,:,irun+1) = data;
    
    %% plot
    figure; setFigParams4Print('landscape')
    colormap gray
    for iplot = 1:irun+1
        subplot(3,2,iplot)
        imagesc(data_mean(:,:,iplot))
        if iplot == irun+1
            title('tun')
        else
            title(num2str(iplot))
        end
    end
    print(fullfile(fnout,[mouse '_' expDate '_check_align']),'-dpdf','-fillpage')
    %% writetiffs
    
    writetiff(data_mean,fullfile(fnout,[mouse '_' expDate '_check_align']))
end
