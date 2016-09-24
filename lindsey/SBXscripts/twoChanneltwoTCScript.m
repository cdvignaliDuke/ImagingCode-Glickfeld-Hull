date = '160624';
ImgFolder = strvcat('005');
time = strvcat('1548');
mouse = 'i916';
nrun = size(ImgFolder,1);

run_str = ['runs-' ImgFolder(1,:)];
if nrun>1
    run_str = [run_str '-' ImgFolder(nrun,:)];
end

data_green = [];
data_red = [];
for irun = 1:nrun
    CD = ['Z:\home\lindsey\Data\2P_images\' date '_' mouse '\' ImgFolder(irun,:)];
    cd(CD);
    imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
    load(imgMatFile);

    nframes = info.config.frames;
    data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
    
    data_red = cat(3, data_red, squeeze(data_temp(2,:,:,:)));
    data_green = cat(3, data_green, squeeze(data_temp(1,:,:,:))); 
end
clear data_temp

data_red_avg = mean(data_red(:,:,1:250),3);
data_green_avg = mean(data_green(:,:,1:250),3);

[out, data_red_reg] = stackRegister(data_red,data_red_avg);
clear data_red
[outs, data_green_reg]=stackRegister_MA(data_green,[],[],out);
clear data_green


data_red_avg = mean(data_red_reg,3);
data_green_avg = mean(data_green_reg,3);

mkdir(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str])) 
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_reg_shifts.mat']), 'out', 'data_red_avg', 'data_green_avg')
clear out outs

for irun = 1:nrun
    fName = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-' mouse '-' date '-' time(irun,:) '.mat'];
    load(fName);
    ds(irun) = input;
end
ds = concatenateDataBlocks(ds);

sz = size(data_green_reg);
if isfield(ds, 'trialOutcomeCell')
    ntrials = length(ds.trialOutcomeCell);
    cFirstStim = celleqel2mat_padded(ds.cFirstStim);
    cTargetOn = celleqel2mat_padded(ds.cTargetOn);
    data_mat = nan(sz(1), sz(2), 100, ntrials);
    for itrial = 1:ntrials
        if ~isnan(cTargetOn(itrial))
            if cTargetOn(itrial) + 29 < sz(3) 
            data_mat(:,:,:,itrial) = data_green_reg(:,:,cFirstStim(itrial)-15:cTargetOn(itrial)+29);
            end
        end
    end
    data_f = nanmean(nanmean(data_mat(:,:,1:15,:),3),4);
    data_df = bsxfun(@minus, data_mat, data_f);
    data_dfof = bsxfun(@rdivide, data_df, data_f);
    dfof_avg = nanmean(data_dfof,4);
    dfof_max_max = mean(dfof_avg(:,:,25:75),3)- mean(dfof_avg(:,:,1:20),3);
else
    ntrials = length(ds.tGratingDirectionDeg);
    nOn = ds.nScansOn;
    nOff = ds.nScansOff;
    data_tr = reshape(data_red_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
    data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
    data_df = bsxfun(@minus, double(data_tr), data_f); 
    data_dfof = bsxfun(@rdivide,data_df, data_f);
    if ds.doTFStim
        TF = celleqel2mat_padded(ds.tGratingTemporalFreqCPS);
        TFs = unique(TF);
        for iTF = 1:length(TFs)
            ind = find(TF == TFs(iTF));
            dfof_max(:,:,iTF) = mean(nanmean(data_dfof(:,:,nOff:nOff+(nOn/2),ind),4),3);
        end
    end  
    dfof_max_max = max(dfof_max,[],3);
end
clear data_mat data_f data_df
for iTF = 1:length(TFs)
	subplot(3,3,iTF)
    imagesc(dfof_max(:,:,iTF))
end
        
bwout = imCellEditInteractive(dfof_max_max);
mask_cell = bwlabel(bwout);


bwout = imCellEditInteractive(data_red_avg);
mask_red = bwlabel(bwout);

mask_np_green = imCellNeuropil(mask_cell, 3, 5);
mask_np_red = imCellNeuropil(mask_red, 3, 5);
save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_mask_cell.mat']), 'mask_red', 'mask_cell', 'mask_np_red', 'mask_np_green', 'dfof_max')

figure; subplot(2,2,1); imagesc(mask_cell)
subplot(2,2,2); imagesc(mask_red)
mask_1 = zeros(size(mask_red));
mask_2 = zeros(size(mask_red));
mask_1(find(mask_red>0)) = 2;
mask_2(find(mask_cell>0)) = 4;
subplot(2,2,3); imagesc(mask_1+mask_2)
subplot(2,2,4); imagesc(bwlabel(mask_1+mask_2))

data_green_reg_down = stackGroupProject(data_green_reg,10);
data_red_tc = stackGetTimeCourses(data_green_reg, mask_red);
data_red_tc_down = stackGetTimeCourses(data_green_reg_down, mask_red);
nCells_red = size(data_red_tc,2);
data_green_tc = stackGetTimeCourses(data_green_reg, mask_cell);
data_green_tc_down = stackGetTimeCourses(data_green_reg_down, mask_cell);
nCells_green = size(data_green_tc,2);

for i = 1:nCells_red
     np_red_tc(:,i) = stackGetTimeCourses(data_green_reg,mask_np_red(:,:,i));
     np_red_tc_down(:,i) = stackGetTimeCourses(data_green_reg_down,mask_np_red(:,:,i));
end
for i = 1:nCells_green
     np_green_tc(:,i) = stackGetTimeCourses(data_green_reg,mask_np_green(:,:,i));
     np_green_tc_down(:,i) = stackGetTimeCourses(data_green_reg_down,mask_np_green(:,:,i));
end

%get weights by maximizing skew
clear x
ii= 0.01:0.01:1;
for i = 1:100
    x(i,:) = skewness(data_red_tc_down-tcRemoveDC(np_red_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_red_tc = data_red_tc-bsxfun(@times,tcRemoveDC(np_red_tc),np_w);

clear x
ii= 0.01:0.01:1;
for i = 1:100
    x(i,:) = skewness(data_green_tc_down-tcRemoveDC(np_green_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_green_tc = data_green_tc-bsxfun(@times,tcRemoveDC(np_green_tc),np_w);

save(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']), 'npSub_green_tc', 'np_green_tc', 'data_green_tc',  'npSub_red_tc', 'np_red_tc', 'data_red_tc')

if isfield(ds, 'trialOutcomeCell')
    data_mat_red = nan(100, nCells_red, ntrials);
    data_mat_green = nan(100, nCells_green, ntrials);
    for itrial = 1:ntrials
        if ~isnan(cTargetOn(itrial))
            if cTargetOn(itrial) + 29 < size(npSub_red_tc,1) 
            data_mat_red(:,:,itrial) = npSub_red_tc(cFirstStim(itrial)-15:cTargetOn(itrial)+29,:);
            data_mat_green(:,:,itrial) = npSub_green_tc(cFirstStim(itrial)-15:cTargetOn(itrial)+29,:);
            end
        end
    end
    data_mat_red = permute(data_mat_red, [2 3 1]);
    data_mat_green = permute(data_mat_green, [2 3 1]);

    red_f = squeeze(nanmean(data_mat_red(:,:,1:15),3));
    red_df = bsxfun(@minus, data_mat_red, red_f);
    red_dfof = bsxfun(@rdivide, red_df, red_f);

    data_red_avg = squeeze(nanmean(red_dfof,2));
    data_red_sem = squeeze(nanstd(red_dfof,[],2))./sqrt(ntrials);

    green_f = squeeze(nanmean(data_mat_green(:,:,1:15),3));
    green_df = bsxfun(@minus, data_mat_green, green_f);
    green_dfof = bsxfun(@rdivide, green_df, green_f);

    data_green_avg = nanmean(green_dfof,2);
    data_green_sem = nanstd(green_dfof,[],2)./sqrt(ntrials);
    
    figure;
    [n n2] = subplotn(nCells_red);
    for i = 1:nCells_red
        subplot(n, n2, i)
        shadedErrorBar(1:100, data_red_avg(i,:), data_red_sem(i,:));
        ylim([-0.05 0.4]) 
        vline([15:11:(11*5)+15])
        vline((11*6)+15,'--k')
        hline(0)
    end
    suptitle('PV-tdTomato Cells')
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_red_cells_avg_all.pdf']), '-dpdf')

    figure;
    [n n2] = subplotn(nCells_green);
    for i = 1:nCells_green
        subplot(n, n2, i)
        shadedErrorBar(1:100, data_green_avg(i,:), data_green_sem(i,:));
        ylim([-0.05 0.4])
        vline([15:11:(11*5)+15])
        vline((11*6)+15,'--k')
        hline(0)
    end
    suptitle('All "responsive" cells')
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_green_cells_avg_all.pdf']), '-dpdf')


    Dir = cell2mat(ds.tGratingDirectionDeg);
    Dirs = unique(Dir);
    for idir = 1:length(Dirs)
        ind = find(Dir == Dirs(idir));
        data_green_dir(:,:,idir) = squeeze(nanmean(green_dfof(:,ind,:),2));
        data_red_dir(:,:,idir) = squeeze(nanmean(red_dfof(:,ind,:),2));
    end

    figure;
    [n n2] = subplotn(nCells_green);
    for i = 1:nCells_green
        for idir = 1:length(Dirs)
            subplot(n, n2, i)
            plot(1:100, data_green_dir(i,:,idir));
            hold on
            ylim([-0.05 0.4])
            vline([15:11:(11*5)+15])
            vline((11*6)+15,'--k')
        end
    end
    legend(num2str(Dirs'));
    suptitle('All "responsive" cells')
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_green_cells_avg_dir.pdf']), '-dpdf')



    figure;
    [n n2] = subplotn(nCells_red);
    for i = 1:nCells_red
        for idir = 1:length(Dirs)
            subplot(n, n2, i)
            plot(1:100, data_red_dir(i,:,idir));
            hold on
            ylim([-0.05 0.4])
            vline([15:11:(11*5)+15])
            vline((11*6)+15,'--k')
        end
    end
    legend(num2str(Dirs'));
    suptitle('PV-tdTomato Cells')
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_red_cells_avg_dir.pdf']), '-dpdf')

else
    nCells_red = size(npSub_red_tc,2);
    data_mat_red = reshape(npSub_red_tc, [nOff+nOn, ntrials, nCells_red]);
    nCells_green = size(npSub_green_tc,2);
    data_mat_green = reshape(npSub_green_tc, [nOff+nOn, ntrials, nCells_green]);
    
    red_f = mean(data_mat_red(nOff/2:nOff,:,:),1);
    red_df = bsxfun(@minus, double(data_mat_red), red_f);
    red_dfof = bsxfun(@rdivide, red_df, red_f);
    
    green_f = mean(data_mat_green(nOff/2:nOff,:,:),1);
    green_df = bsxfun(@minus, double(data_mat_green), green_f);
    green_dfof = bsxfun(@rdivide, green_df, green_f);
    
    data_red_avg = squeeze(nanmean(red_dfof,2));
    data_red_sem = squeeze(nanstd(red_dfof,[],2))./sqrt(ntrials);
    
    data_green_avg = squeeze(nanmean(green_dfof,2));
    data_green_sem = squeeze(nanstd(green_dfof,[],2)./sqrt(ntrials));
    
    figure;
    [n n2] = subplotn(nCells_red);
    for i = 1:nCells_red
        subplot(n, n2, i)
        shadedErrorBar(1:size(red_dfof,1), data_red_avg(:,i), data_red_sem(:,i));
        ylim([-0.05 0.05]) 
        vline(nOff+1)
        hline(0)
    end
    suptitle('PV-tdTomato Cells')
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_red_cells_avg_all.pdf']), '-dpdf')

    figure;
    [n n2] = subplotn(nCells_green);
    for i = 1:nCells_green
        subplot(n, n2, i)
        shadedErrorBar(1:size(green_dfof,1), data_green_avg(:,i), data_green_sem(:,i));
        ylim([-0.05 0.2]) 
        vline(nOff+1)
        hline(0)
    end
    suptitle('All "responsive" cells')
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_green_cells_avg_all.pdf']), '-dpdf')
    
    
    for iTF = 1:length(TFs)
        ind = find(TF == TFs(iTF));
        data_green_TF(:,:,iTF) = squeeze(nanmean(green_dfof(:,ind,:),2));
        data_red_TF(:,:,iTF) = squeeze(nanmean(red_dfof(:,ind,:),2));
    end

    figure;
    [n n2] = subplotn(nCells_green);
    for i = 1:nCells_green
        for iTF = 1:length(TFs)
            subplot(n, n2, i)
            plot(1:size(green_dfof,1), data_green_TF(:,i,iTF));
            hold on
            ylim([-0.1 0.6]) 
            vline(nOff+1)
            hline(0)
        end
    end
    suptitle('All "responsive" cells')
    legend(num2str(TFs'));
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_green_cells_avg_TF.pdf']), '-dpdf')



   figure;
    [n n2] = subplotn(nCells_red);
    for i = 1:nCells_red
        for iTF = 1:length(TFs)
            subplot(n, n2, i)
            plot(1:size(red_dfof,1), data_red_TF(:,i,iTF));
            hold on
            ylim([-0.05 0.1]) 
            vline(nOff+1)
            hline(0)
        end
    end
    suptitle('PV-tdTomato Cells')
    legend(num2str(TFs'));
    print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_red_cells_avg_TF.pdf']), '-dpdf')
end








