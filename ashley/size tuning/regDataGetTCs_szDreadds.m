clear all
close all
ds = 'szTuning_dreadds_AL';
rc = behavConstsAV;
eval(ds)
slct_expt = 1:2;
%%
for iexp = slct_expt

    mouse = expt(iexp).mouse;
    subnum = mouse;
    expDate = expt(iexp).date;
    runs = expt(iexp).sizeTuningFolder;
    nrun = length(runs);

    fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
    mkdir(fnout)
    %% register all size tuning experiments
    for irun = 1:nrun
        dataFolder = runs{irun};
        fName = [dataFolder '_000_000'];
        switch expt(iexp).saveLoc
            case 'ashley'
                data = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName);
            case 'kevin'
                fdir = ['Y:\home\kevin\Data\2P\' expDate '_i' mouse '\' dataFolder];
                data = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName,[],fdir);
            otherwise
                error('identify data source')
        end
        if irun == 1
            nfr = size(data,3); 

            nMeanImages = 9;
            [nRows,nCols] = optimizeSubplotDim(nMeanImages+1);

            randStarterFrames = sort(randsample(nfr-100,nMeanImages));
            setFigParams4Print('landscape')
            figure
            suptitle([subnum '-' expDate])
            for iimg = 1:nMeanImages
                subplot(nRows,nCols,iimg)
                imagesc(mean(data(:,:,randStarterFrames(iimg):randStarterFrames(iimg)+100),3))
                title(sprintf('%s:%s',num2str(randStarterFrames(iimg)),num2str(randStarterFrames(iimg)+100)))
            end
            print(['Z:\Analysis\Expt Summaries\' ds '\rand image samples_' subnum '-' expDate],'-dpdf','-fillpage')
            savefig(['Z:\Analysis\Expt Summaries\' ds '\rand image samples_' subnum '-' expDate])

            regImageFramesStart = input('Enter frame number for registration:');
            regImage = mean(data(:,:,regImageFramesStart:(regImageFramesStart+99)),3);
            save(fullfile(fnout,'regImage'),'regImage')
        else
            load(fullfile(fnout,'regImage'))
        end

        [outs,data_reg] = stackRegister(data,regImage);
        mkdir(fullfile(fnout,dataFolder))
        save(fullfile(fnout,dataFolder,'registration'),'outs','regImage')
        clear data

        if irun == 1
            F0 = uint16(mean(data_reg,3));
            dFF = (data_reg - F0)./F0;
            maxDFF = max(dFF,[],3);
            writetiff(maxDFF,fullfile(fnout,'maxDFF_cellSelectionImage'))

            bw = imCellEditInteractive(maxDFF);
            green_mask = bwlabel(bw);
            save(fullfile(fnout,'cellMask'),'green_mask');   
        else
            load(fullfile(fnout,'cellMask'))     
        end

        tc = stackGetTimeCourses(data_reg,green_mask);
        dataTC = getWeightedNeuropilTimeCourse(data_reg,tc,green_mask,4,6);
        save(fullfile(fnout,dataFolder,'timecourses'),'dataTC')
        clear data_reg tc dataTC
    end

    %% register red channel and save image
    if ~isempty(expt(iexp).labelFolder)
        labelFolder = expt(iexp).labelFolder;
        fName = [labelFolder '_000_000'];
        switch expt(iexp).saveLoc
            case 'ashley'
                data = loadsbx_choosepmt(2,mouse,expDate,labelFolder,fName);
            case 'kevin'
                fdir = ['Y:\home\kevin\Data\2P\' expDate '_i' mouse '\' labelFolder];
                data = loadsbx_choosepmt(2,mouse,expDate,labelFolder,fName,[],fdir);
            otherwise
                error('identify data source')
        end
        figure;imagesc(mean(data,3))
        [outs,data_red_reg] = stackRegister(data,regImage);
        data_red_reg_down = stackGroupProject(data_red_reg,10);
        figure;imagesc(mean(data_red_reg_down,3))
        writetiff(data_red_reg_down,fullfile(fnout,'mCherryImage'))
        clear data data_red_reg data_red_reg_down
    end
    %% register retinotopy 
    load(fullfile(fnout,'cellMask'))
    load(fullfile(fnout,'regImage'))
    dataFolder = expt(iexp).retinotopyFolder{1};
    fName = [dataFolder '_000_000'];
    switch expt(iexp).saveLoc
        case 'ashley'
            data = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName);
        case 'kevin'
            fdir = ['Y:\home\kevin\Data\2P\' expDate '_i' mouse '\' dataFolder];
            data = loadsbx_choosepmt(1,mouse,expDate,dataFolder,fName,[],fdir);
        otherwise
            error('identify data source')
    end
    [outs,data_reg] = stackRegister(data,regImage);
    mkdir(fullfile(fnout,dataFolder))
    save(fullfile(fnout,dataFolder,'registration'),'outs','regImage')
    clear data
    tc = stackGetTimeCourses(data_reg,green_mask);
    dataTC = getWeightedNeuropilTimeCourse(data_reg,tc,green_mask,4,6);
    save(fullfile(fnout,dataFolder,'timecourses'),'dataTC')
    clear data_reg tc dataTC
end