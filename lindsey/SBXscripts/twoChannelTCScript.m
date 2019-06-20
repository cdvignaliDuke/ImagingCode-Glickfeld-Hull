date = '190214';
ImgFolder = strvcat('005');
RedFolder = strvcat('006');
time = strvcat('1612');
mouse = '1211';
nrun = size(ImgFolder,1);
irun = 1;
fout = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P', [date '_i' mouse], [date '_i' mouse '_' ImgFolder]);
if ~exist(fout)
    mkdir(fout)
end

run_str = ['runs-' ImgFolder(1,:)];
if nrun>1
    run_str = [run_str '-' ImgFolder(nrun,:)];
end
clear global
CD = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\ashley\data\' mouse '\two-photon imaging\' date '\' ImgFolder(irun,:)];
cd(CD);
imgMatFile = [ImgFolder(irun,:) '_000_000.mat'];
load(imgMatFile);

nframes = info.config.frames;
data_temp = sbxread([ImgFolder(irun,:) '_000_000'],0,nframes);
data_green = squeeze(data_temp(1,:,:,:)); 
clear data_temp

data_avg = mean(data_green(:,:,5001:5500),3);

[out, data_reg] = stackRegister(data_green,data_avg);
save(fullfile(fout, [date '_' mouse '_' run_str '_reg.mat']), 'out','data_avg')
clear data_green

clear global
CD = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\ashley\data\' mouse '\two-photon imaging\' date '\' RedFolder(irun,:)];
cd(CD);
imgMatFile = [RedFolder(irun,:) '_000_000.mat'];
load(imgMatFile);

nframes = info.config.frames;
data_temp = sbxread([RedFolder(irun,:) '_000_000'],0,nframes);
ref_green = mean(squeeze(data_temp(1,:,:,:)),3);
ref_red = mean(squeeze(data_temp(2,:,:,:)),3); 
clear data_temp

[out, data_g] = stackRegister(ref_green,data_avg);
[outs, data_r] = stackRegister_MA(ref_red,[],[],out);

%% segment cells
fName = ['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data\data-i' mouse '-' date '-' time(irun,:) '.mat'];
load(fName);


ntrials = length(input.tGratingDirectionDeg);
nOn = input.nScansOn;
nOff = input.nScansOff;
nframes = (nOn+nOff) .* ntrials;
data_reg = data_reg(:,:,1:nframes);

sz = size(data_reg);
data_tr = reshape(data_reg,[sz(1), sz(2), nOn+nOff, ntrials]);
data_f = mean(data_tr(:,:,nOff/2:nOff,:),3);
data_df = bsxfun(@minus, double(data_tr), data_f); 
data_dfof = bsxfun(@rdivide,data_df, data_f);
con_mat = celleqel2mat_padded(input.tGratingContrast);
cons = unique(con_mat);
ncon = length(cons);
dfof_all = zeros(sz(1),sz(2),ncon);
for icon = 1:ncon
    ind = find(con_mat == cons(icon));
    dfof_all(:,:,icon) = mean(nanmean(data_dfof(:,:,nOff:nOff+nOn,ind),4),3);
end
clear data_tr data_f data_df data_dfof

img_all = cat(3,data_r,data_r,dfof_all);
    
mask_exp = zeros(sz(1),sz(2));
mask_all = zeros(sz(1), sz(2));
mask_data = img_all;

for iStim = 1:size(img_all,3)
    mask_data_temp = mask_data(:,:,iStim);
    mask_data_temp(find(mask_exp >= 1)) = 0;
    bwout = imCellEditInteractive(mask_data_temp);
    mask_all = mask_all+bwout;
    if iStim == 2
        mask_red = mask_all;
    end
    mask_exp = imCellBuffer(mask_all,3)+mask_all;
    close all
end
mask_cell= bwlabel(mask_all);
figure; imagesc(mask_cell)

red_cells = unique(mask_cell(find(mask_red)));

mask_np = imCellNeuropil(mask_cell, 3, 5);
save(fullfile(fout, [date '_' mouse '_' run_str '_mask_cell.mat']), 'img_all', 'mask_cell', 'mask_np', 'red_cells')


%% neuropil subtraction
down = 5;

data_tc = stackGetTimeCourses(data_reg, mask_cell);
data_reg_down  = stackGroupProject(data_reg,down);
data_tc_down = stackGetTimeCourses(data_reg_down, mask_cell);
nCells = size(data_tc,2);
np_tc = zeros(sz(3),nCells);
np_tc_down = zeros(floor(sz(3)./down), nCells);
for i = 1:nCells
     np_tc(:,i) = stackGetTimeCourses(data_reg,mask_np(:,:,i));
     np_tc_down(:,i) = stackGetTimeCourses(data_reg_down,mask_np(:,:,i));
     fprintf(['Cell #' num2str(i) '%s/n']) 
end
%get weights by maximizing skew
ii= 0.01:0.01:1;
x = zeros(length(ii), nCells);
for i = 1:100
    x(i,:) = skewness(data_tc_down-tcRemoveDC(np_tc_down*ii(i)));
end
[max_skew ind] =  max(x,[],1);
np_w = 0.01*ind;
npSub_tc = data_tc-bsxfun(@times,tcRemoveDC(np_tc),np_w);
clear data_reg data_reg_down

save(fullfile(fout, [date '_' mouse '_' run_str '_TCs.mat']), 'data_tc', 'np_tc', 'npSub_tc')

%%
red_int = zeros(1,nCells);
green_int = zeros(1,nCells);
red_np = zeros(1,nCells);
for iCell = 1:nCells
    ind = find(mask_cell == iCell);
    red_int(1,iCell) = mean(data_r(ind));
    green_int(1,iCell) = mean(data_avg(ind));
    ind_np = find(mask_np(:,:,iCell));
    red_np(1,iCell) = mean(data_r(ind_np));
end
figure; 
subplot(2,2,1)
hist(red_int(red_cells)./green_int(red_cells))
title('Red/Green- SOM+')
xlim([0 2])
subplot(2,2,2)
hold on
hist(red_int(setdiff(1:nCells,red_cells))./green_int(setdiff(1:nCells,red_cells)))
title('Red/Green- SOM-')
xlim([0 2])  
subplot(2,2,3)
hist(red_np(red_cells)./red_int(red_cells),[0:0.1:2])
title('NP/Int- SOM+')
xlim([0 2]) 
ylim([0 15])
subplot(2,2,4)
hist(red_np(setdiff(1:nCells,red_cells))./red_int(setdiff(1:nCells,red_cells)),[0:0.1:2])
title('NP/Int- SOM-')
xlim([0 2])
ylim([0 15])
print(fullfile(fout, [date '_i' mouse '_' ImgFolder '_intensityHists.pdf']),'-dpdf','-bestfit');


%% contrast response
%trial_tc = permute(reshape(npSub_tc, [nOn+nOff ntrials nCells]),[1 3 2]);
trial_tc = nan(nOff+nOn+nOff, nCells, ntrials);
for itrial = 1:ntrials
    if ((itrial-1).*(nOn+nOff))+nOff+nOn+nOff <= nframes
        trial_tc(:,:,itrial) = npSub_tc(((itrial-1).*(nOn+nOff))+1:((itrial-1).*(nOn+nOff))+nOff+nOn+nOff,:);
    end
end
tc_f = mean(trial_tc(0.75.*nOff:nOff,:,:),1);
tc_dfof = (trial_tc-tc_f)./tc_f;
dfof_con = zeros(nOff+nOn+nOff, nCells, ncon);
figure;
[n n2] = subplotn(nCells);
for icon = 1:ncon
    ind = find(con_mat == cons(icon));
    dfof_con(:,:,icon) = nanmean(tc_dfof(:,:,ind),3);
    for iCell = 1:nCells
        subplot(n,n2,iCell)
        plot(dfof_con(:,iCell,icon))
        hold on
        if find(red_cells == iCell)
            title('SOM')
        end
    end
end
suptitle([mouse ' ' date ' Contrast response'])
print(fullfile(fout, [date '_i' mouse '_' ImgFolder '_contrastTCs.pdf']),'-dpdf','-bestfit');

con_resp_avg = zeros(nCells, ncon, 2);
[n n2] = subplotn(nCells);
for icon = 1:ncon
    ind = find(con_mat == cons(icon));
    con_resp_avg(:,icon,1) = squeeze(nanmean(mean(tc_dfof(nOff+5:nOff+nOn+5,:,ind),1),3))';
    con_resp_avg(:,icon,2) = squeeze(nanstd(mean(tc_dfof(nOff+5:nOff+nOn+5,:,ind),1),[],3)./sqrt(length(ind)))';
end
figure;
for iCell = 1:nCells
    subplot(n,n2,iCell)
    errorbar(cons,con_resp_avg(iCell,:,1),con_resp_avg(iCell,:,2),'-o')
    hold on
    if find(red_cells == iCell)
        title('SOM')
    end
end
suptitle([mouse ' ' date ' Contrast response avg'])
print(fullfile(fout, [date '_i' mouse '_' ImgFolder '_contrastResp.pdf']),'-dpdf','-bestfit');

figure;
subplot(1,2,1)
errorbar(cons, mean(con_resp_avg(setdiff(1:nCells,red_cells),:,1),1)', std(con_resp_avg(setdiff(1:nCells,red_cells),:,1),[],1)./sqrt(length(setdiff(1:nCells,red_cells)))','-o')
hold on
errorbar(cons, mean(con_resp_avg(red_cells,:,1),1)', std(con_resp_avg(red_cells,:,1),[],1)./sqrt(length(red_cells))','-o')
subplot(1,2,2)
con_resp_norm = con_resp_avg(:,:,1)./max(con_resp_avg(:,:,1),[],2);
errorbar(cons, mean(con_resp_norm(setdiff(1:nCells,red_cells),:),1), std(con_resp_norm(setdiff(1:nCells,red_cells),:),[],1)./sqrt(length(setdiff(1:nCells,red_cells))),'-o')
hold on
errorbar(cons, mean(con_resp_norm(red_cells,:),1), std(con_resp_norm(red_cells,:),[],1)./sqrt(length(red_cells)),'-o')

legend({'SOM-', 'SOM+'},'location','southeast')

%remove the highest con (error in code to 1.6 instead of 1)
figure;
subplot(1,2,1)
errorbar(cons(1:5).*100, mean(con_resp_avg(setdiff(1:nCells,red_cells),1:5,1),1)', std(con_resp_avg(setdiff(1:nCells,red_cells),1:5,1),[],1)./sqrt(length(setdiff(1:nCells,red_cells)))','-o')
hold on
errorbar(cons(1:5).*100, mean(con_resp_avg(red_cells,1:5,1),1)', std(con_resp_avg(red_cells,1:5,1),[],1)./sqrt(length(red_cells))','-o')
xlim([0 100])
ylabel('dF/F')
xlabel ('Contrast')
subplot(1,2,2)
con_resp_norm = con_resp_avg(:,1:5,1)./max(con_resp_avg(:,1:5,1),[],2);
errorbar(cons(1:5).*100, mean(con_resp_norm(setdiff(1:nCells,red_cells),:),1), std(con_resp_norm(setdiff(1:nCells,red_cells),:),[],1)./sqrt(length(setdiff(1:nCells,red_cells))),'-o')
hold on
errorbar(cons(1:5).*100, mean(con_resp_norm(red_cells,:),1), std(con_resp_norm(red_cells,:),[],1)./sqrt(length(red_cells)),'-o')
xlim([0 100])
ylabel('Norm dF/F')
xlabel ('Contrast')
legend({'SOM-', 'SOM+'},'location','southeast')
print(fullfile(fout, [date '_i' mouse '_' ImgFolder '_avgContrastResp.pdf']),'-dpdf','-bestfit');

%% contrast fits
conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));
opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
con_resp_all = [];
con_list = [];
for icon = 1:ncon-1
    ind = find(con_mat == cons(icon));
    con_resp_all = [con_resp_all; squeeze(mean(tc_dfof(nOff+5:nOff+nOn+5,:,ind),1))'];
    con_list = [con_list; cons(icon).*ones(length(ind),1)];
end
ind = find(isnan(con_resp_all(:,1)));
con_resp_all(ind,:) = [];
con_list(ind,:) = [];

lb = [0 0 0.1 1];
ub = [Inf Inf 0.8 Inf];
cF = zeros(4,nCells);
res = zeros(1,nCells);
conRng = 0:0.01:1;
figure;
for iCell = 1:nCells
    cRi = con_resp_all(:,iCell);
    SStot = sum((cRi-mean(cRi)).^2);
    [min_val, min_ind] = min(abs(0.5- (con_resp_avg(iCell,1:5,1)./max(con_resp_avg(iCell,1:5,1)))),[],2);
    x0 = [0 mean(cRi) cons(min_ind) 1]; %BL Rmax C50 n
    [cF(:,iCell), res(1,iCell)] = lsqcurvefit(conModelH,x0,con_list,cRi,lb,ub,opts);
    subplot(n,n2,iCell)
    plot(cons(1:5),con_resp_avg(iCell,1:5,1),'o')
    hold on
    if find(red_cells == iCell)
        title('SOM')
    end
    fitout = conModelH(cF(:,iCell),conRng);
    plot(conRng,fitout,'-r');
end



rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = data_r./max(data_r(:));
rgb(:,:,2) = data_avg./max(data_avg(:));

figure;
imagesc(rgb(200:200+sz(1)/2,:,:))
truesize
hold on
for iC = 1:nCells
    [i,j] = find(mask_cell(200:200+sz(1)/2,:) == iC);
    c = [mean(i,1) mean(j,1)];
    text(c(2),c(1),num2str(iC))
end

ex_cells = [4 29 9 20];
figure; 
imagesc(rgb(200:200+sz(1)/2,:,:))
truesize
hold on
[i1,j1] = find(mask_cell(200:200+sz(1)/2,:) == ex_cells(1));
[i2,j2] = find(mask_cell(200:200+sz(1)/2,:) == ex_cells(2));
[i3,j3] = find(mask_cell(200:200+sz(1)/2,:) == ex_cells(3));
[i4,j4] = find(mask_cell(200:200+sz(1)/2,:) == ex_cells(4));
c1 = [mean(i1,1) mean(j1,1)];
c2 = [mean(i2,1) mean(j2,1)];
c3 = [mean(i3,1) mean(j3,1)];
c4 = [mean(i4,1) mean(j4,1)];
scatter(c1(2),c1(1),'og')
scatter(c2(2),c2(1),'og')
scatter(c3(2),c3(1),'or')
scatter(c4(2),c4(1),'or')
print(fullfile(fout, [date '_i' mouse '_' ImgFolder '_RGB_FOV_exCells.pdf']),'-dpdf','-bestfit');

figure;
for iCell = 1:4
    iC = ex_cells(iCell);
    subplot(2,2,iCell)
    errorbar(cons(1:5),con_resp_avg(iC,1:5,1),con_resp_avg(iC,1:5,2),'o')
    hold on
    fitout = conModelH(cF(:,iC),conRng);
    plot(conRng,fitout,'-r');
    ylim([-0.02 0.15])
    if find(red_cells == iC)
        title('SOM')
    end
    ylabel('dF/F')
    xlabel('Contrast')
end
print(fullfile(fout, [date '_i' mouse '_' ImgFolder '_conFits_exCells.pdf']),'-dpdf','-bestfit');




%%

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








