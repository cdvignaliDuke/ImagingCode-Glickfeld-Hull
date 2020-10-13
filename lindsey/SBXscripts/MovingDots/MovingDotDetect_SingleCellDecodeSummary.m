close all;clear all;clc;

ds = 'MovingDotDetect_ExptList';
rc = behavConstsAV;
eval(ds)

area_list = {'V1', 'AL', 'PM'};
singleCellPctCorrect = cell(size(area_list));
resp_cells = cell(size(area_list));
inc_resp_cells = cell(size(area_list));
tuned_cells = cell(size(area_list));
fract_resp = cell(size(area_list));
fract_inc_resp = cell(size(area_list));
fract_tuned = cell(size(area_list));
max_spd = cell(size(area_list));
max_spd_nobase = cell(size(area_list));
max_resp = cell(size(area_list));
cell_weights = cell(size(area_list));
max_spd_grating = cell(size(area_list));
TF_resp_cells = cell(size(area_list));
TF_tuned_cells = cell(size(area_list));
dfof_spd_resp_all = cell(size(area_list));
spd_resp_all = cell(size(area_list));
expt_id = cell(size(area_list));
mice = cell(size(expt));

for iexp = 1:length(expt)
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).dotFolder;
    nrun = length(ImgFolder);
    dot_run_str = catRunName(cell2mat(ImgFolder), nrun);
    ImgFolder = expt(iexp).gratFolder;
    nrun = length(ImgFolder);
    grating_run_str = catRunName(cell2mat(ImgFolder), nrun);
    fprintf(['Expt ' num2str(iexp) ': ' date ' ' mouse '\n'])
    mice{iexp} = [expt(iexp).mouse '_'];
    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
    
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dot_run_str], [date '_' mouse '_' dot_run_str '_cellResp.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dot_run_str], [date '_' mouse '_' dot_run_str '_trResp.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dot_run_str], [date '_' mouse '_' dot_run_str '_singleCellDecode.mat']))
    
    iarea = find(strcmp(area_list, area));
    resp_cells{iarea} = [resp_cells{iarea}; resp_ind+length(singleCellPctCorrect{iarea})];
    inc_resp_cells{iarea} = [inc_resp_cells{iarea};  inc_resp_ind+length(singleCellPctCorrect{iarea})];
    tuned_cells{iarea} = [tuned_cells{iarea} anova_ind+length(singleCellPctCorrect{iarea})];
    fract_resp{iarea} = [fract_resp{iarea} length(resp_ind)./length(pctCorrectTarget_ho_singlecell)];
    fract_inc_resp{iarea} = [fract_inc_resp{iarea} length(inc_resp_ind)./length(resp_ind)];
    fract_tuned{iarea} = [fract_tuned{iarea} length(anova_ind)./length(inc_resp_ind)];
    singleCellPctCorrect{iarea} = [singleCellPctCorrect{iarea} pctCorrectTarget_ho_singlecell];
    max_spd{iarea} = [max_spd{iarea}; max_ind];
    max_spd{iarea}(find(max_val<0)) = 0;
    max_resp{iarea} = [max_resp{iarea}; max_val];
    cell_weights{iarea} = [cell_weights{iarea} cellW];
    max_spd_nobase{iarea} = [max_spd_nobase{iarea}; max_ind_nobase];
    expt_id{iarea} = [expt_id{iarea} iexp.*ones(size(pctCorrectTarget_ho_singlecell))];
    
    dfof_spd_resp_all{iarea} = cat(2, dfof_spd_resp_all{iarea}, data_dfof_spd);
    
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' grating_run_str], [date '_' mouse '_' grating_run_str '_dfofData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' grating_run_str], [date '_' mouse '_' grating_run_str '_stimData.mat']))
    max_spd_grating{iarea} = max_ind_tf;
    TF_resp_cells{iarea} = TFresp_ind;
    TF_tuned_cells{iarea} = TFtuned_ind;
    
    spd_resp_all{iarea} = cat(1, spd_resp_all{iarea}, spd_resp_mat);
end
mouse_list = cell2mat(unique(mice));
summary_pn = fullfile(LG_base, '\Analysis\2P\SpeedTuning');
%% responsive cells
% for iarea=2:length(area_list)
%     fprintf([area_list{iarea} ':\n Responsive: ' num2str(length(resp_cells{iarea})) '/' num2str(length(singleCellPctCorrect{iarea})) '- ' ...
%         num2str(mean(fract_resp{iarea},2)) '+/-' num2str(num2str(std(fract_resp{iarea},[],2)./sqrt(length(fract_resp{iarea}))))...
%         '\n Inc Responsive: ' num2str(length(inc_resp_cells{iarea})) '/' num2str(length(resp_cells{iarea})) '- ' ...
%         num2str(mean(fract_inc_resp{iarea},2)) '+/-' num2str(num2str(std(fract_inc_resp{iarea},[],2)./sqrt(length(fract_inc_resp{iarea}))))...
%         '\n Tuned: ' num2str(length(tuned_cells{iarea})) '/' num2str(length(inc_resp_cells{iarea})) '- ' ...
%         num2str(mean(fract_tuned{iarea},2)) '+/-' num2str(num2str(std(fract_tuned{iarea},[],2)./sqrt(length(fract_tuned{iarea})))) '\n']);
% end
 
resp=[size(dfof_spd_resp_all{2},2) size(resp_cells{2},1); size(dfof_spd_resp_all{3},2) size(resp_cells{3},1)];
 e = sum(resp,2)*sum(resp)/sum(resp(:)); % expected
 X2 = (resp-e).^2./e;
 X2 = sum(X2(:)); % chi square
 df = prod(size(resp)-[1 1]); % degree of freedom
 p_resp = 1-chi2cdf(X2,df);
 [AL_m AL_CI] = binofit(size(resp_cells{2},1), size(dfof_spd_resp_all{2},2));
 [PM_m PM_CI] = binofit(size(resp_cells{3},1), size(dfof_spd_resp_all{3},2));
 fprintf(['Responsive CI: AL- ' num2str(chop(AL_CI,2)) '; PM- ' num2str(chop(PM_CI,2)) '\n' ]);
fprintf(['Responsive: ' num2str(chop(p_resp,2)) '\n']);

inc=[size(resp_cells{2},1) size(inc_resp_cells{2},1); size(resp_cells{3},1) size(inc_resp_cells{3},1)];
 e = sum(inc,2)*sum(inc)/sum(inc(:)); % expected
 X2 = (inc-e).^2./e;
 X2 = sum(X2(:)); % chi square
 df = prod(size(inc)-[1 1]); % degree of freedom
 p_inc = 1-chi2cdf(X2,df);
 [AL_m AL_CI] = binofit(size(inc_resp_cells{2},1), size(resp_cells{2},1));
 [PM_m PM_CI] = binofit(size(inc_resp_cells{3},1), size(resp_cells{3},1));
 fprintf(['Inc Responsive CI: AL- ' num2str(chop(AL_CI,2)) '; PM- ' num2str(chop(PM_CI,2)) '\n' ]);
fprintf(['Inc Responsive: ' num2str(chop(p_inc,2)) '\n']);

tuned=[size(inc_resp_cells{2},1) size(tuned_cells{2},2);size(inc_resp_cells{3},1) size(tuned_cells{3},2)];
 e = sum(tuned,2)*sum(tuned)/sum(tuned(:)); % expected
 X2 = (tuned-e).^2./e;
 X2 = sum(X2(:)); % chi square
 df = prod(size(tuned)-[1 1]); % degree of freedom
 p_tuned = 1-chi2cdf(X2,df);

 [AL_m AL_CI] = binofit(size(tuned_cells{2},2), size(inc_resp_cells{2},1));
 [PM_m PM_CI] = binofit(size(tuned_cells{3},2), size(inc_resp_cells{3},1));
 fprintf(['Tuned CI: AL- ' num2str(chop(AL_CI,2)) '; PM- ' num2str(chop(PM_CI,2)) '\n' ]);
fprintf(['Tuned: ' num2str(chop(p_tuned,2)) '\n']);
%% average response to each speed
nSpd = size(dfof_spd_resp_all{1},4);
tt = [-22:97].*(1000/frameRateHz);
figure;
for iarea=2:length(area_list)
    for iSpd = 1:nSpd
        subplot(2,nSpd,iSpd+((iarea-2).*nSpd))
        shadedErrorBar(tt, mean(dfof_spd_resp_all{iarea}(:,inc_resp_cells{iarea},2,iSpd),2), std(dfof_spd_resp_all{iarea}(:,inc_resp_cells{iarea},2,iSpd),[],2)./sqrt(length(inc_resp_cells{iarea})));
        ylim([-0.01 .2])
        xlim([-500 3000])
        vline([200 733], 'k')
        if iarea == 2
            title([num2str(spd_step(iSpd+1)) ' deg/s'])
        end
        if iSpd == 1
            ylabel(area_list{iarea})
        end
    end
end
print(fullfile(summary_pn, [mouse_list 'speedSummary_avgTCs.pdf']),'-dpdf','-fillpage')
%% fit speed response curves
spdRng = 0.5:0.01:30;
nSpd = length(spd_step(2:end));
sn = zeros(1,nSpd);
s = zeros(1);
conModelH = @(coefs,cdata) coefs(1) + coefs(2)*(cdata.^coefs(4))./(cdata.^coefs(4)+coefs(3).^coefs(4));
opts = optimoptions('lsqcurvefit','Display','off'); %,'Algorithm','levenberg-marquardt'
for iarea = 1:length(area_list)
    fprintf([area_list{iarea} '\n'])
    spdSummary(iarea).spdFit = struct('resp',sn,'fit',sn,'C50r',s,'Rsq',s,'x0',sn);
    nCells = size(spd_resp_all{iarea},1);
    nCon = size(spd_resp_all{iarea},2);
    spdSummary(iarea).spdFit(nCells,nCon) = spdSummary(iarea).spdFit;
    for iCell = 1:nCells
        for iCon = 1:nCon
            fprintf(['BaseCon = ' num2str(iCon-1) '\n'])
            cRi = squeeze(spd_resp_all{iarea}(iCell,iCon,:))';
            lb = [0 0 0.1*spd_step(end) 1];
            ub = [Inf Inf 0.8*spd_step(end) Inf];
            SStot = sum((cRi-mean(cRi)).^2);
            R2best = -Inf;
            for i=1%1:4
                x0 = [cRi(1) max(cRi) 0.5*spd_step(end) 3]; %BL Rmax C50 n
                [cF, res] = lsqcurvefit(conModelH,x0,spd_step(2:end),cRi,lb,ub,opts);
                R2 = 1-res/SStot;
                if R2>R2best
                    R2best = R2;
                    cFbest = cF;
                    x0best = x0;
                end
            end
            cF = cFbest;
            R2 = R2best;
            spdSummary(iarea).spdFit(iCell,iCon).fit = cF;
            spdSummary(iarea).spdFit(iCell,iCon).Rsq = R2;
            spdSummary(iarea).spdFit(iCell,iCon).x0 = x0best;

            fitout = conModelH(cF,spdRng);
            R50 = fitout(1)+(fitout(end)-fitout(1))/2;
            fitout50rect = abs(fitout - R50);
            i50 = find(fitout50rect == min(fitout50rect),1);
            C50 = spdRng(i50);
            spdSummary(iarea).spdFit(iCell,iCon).C50r = C50;

            fprintf('Cell %d/%d fit: BL=%.3f Rmax=%.3f C50=%.3f n=%.2f : Rsq=%.3f C50r=%.3f\n',iCell,nCells,cF(1),cF(2),cF(3),cF(4),R2,C50)

        end
    end
end
fprintf('Done, saving...\n')
save(fullfile(summary_pn, [mouse_list 'speedFitSummary.mat']),'spdSummary');
%%
load(fullfile(summary_pn, [mouse_list 'speedFitSummary.mat']))
fract_wellfit = cell(size(area_list));
for iarea = 1:length(area_list)
    expts = unique(expt_id{iarea});
    fract_wellfit{iarea} = zeros(1,length(expts));
    R2 = [spdSummary(iarea).spdFit(:,2).Rsq];
    fprintf([area_list{iarea} '- R2>0.5: ' num2str(length(find(R2(tuned_cells{iarea})>0.5))) '/' num2str(length(tuned_cells{iarea})) '\n'])
    for iexpt = 1:length(expts)
        ind = intersect(tuned_cells{iarea},find(expt_id{iarea} == expts(iexpt)));
        fract_wellfit{iarea}(1,iexpt) = length(find(R2(ind)>0.5))./length(ind);
    end
    fprintf([num2str(mean(fract_wellfit{iarea},2)) '+/-' num2str(num2str(std(fract_wellfit{iarea},[],2)./sqrt(length(fract_wellfit{iarea})))) '\n']);
end    

%%

figure;
start = 0;
C50_sub = cell(size(area_list));
C50_sub_all = [];
area_all = [];
n_sub_all = [];
for iarea = 1:length(area_list)
    subplot(4,3,start+1)
    R2 = [spdSummary(iarea).spdFit(tuned_cells{iarea},2).Rsq];
    R2(find(R2<0))=0;
    edges = [0:0.1:1];
    hist(R2,edges)
    xlabel('Rsq')
    ylabel('# Cells')
    title(area_list{iarea})
    subplot(4,3,(start+2))
    C50 = [spdSummary(iarea).spdFit(tuned_cells{iarea},2).C50r];
    C50_sub{iarea} = C50(find(R2>0.95));
    edges = [0:1:30];
    hist(C50_sub{iarea},edges)
    xlabel('C50- R2>0.5')
    ylabel('# Cells')
    title([num2str(length(find(R2>0.5))) '/' num2str(length(R2)) '- ' num2str(chop(length(find(R2>0.5))./length(R2),2))])
    subplot(4,3,(start+3))
    n = reshape([spdSummary(iarea).spdFit(tuned_cells{iarea},2).fit],[4, length(tuned_cells{iarea})]);
    n_sub{iarea} = n(4,find(R2>0.5));
    edges = [0:10:200];
    hist(n_sub{iarea},edges)
    xlabel('exponent- R2>0.5')
    ylabel('# Cells')
    subplot(4,2,7)
    errorbar(iarea, median(C50_sub{iarea},2), std(C50_sub{iarea},[],2),'o')
    hold on
    C50_sub_all = [C50_sub_all; C50_sub{iarea}'];
    area_all = [area_all; iarea.*ones(size(C50_sub{iarea}))'];
    ylabel('C50')
    set(gca,'XTick', 1:3, 'XTickLabel',area_list)
    xlim([0 4])
    ylim([0 25])
    subplot(4,2,8)
    errorbar(iarea, median(n_sub{iarea},2), std(n_sub{iarea},[],2),'o')
    hold on
    n_sub_all = [n_sub_all; n_sub{iarea}'];
    ylabel('Exponent')
    set(gca,'XTick', 1:3, 'XTickLabel',area_list)
    xlim([0 4])
    ylim([0 100])
    start = start+3;
end

print(fullfile(summary_pn, [mouse_list 'speedSummary_Fits.pdf']),'-dpdf','-fillpage')

figure;
iarea = 1;
R2 = [spdSummary(iarea).spdFit(tuned_cells{iarea},2).Rsq];
R2(find(R2<0))=0;
ind = find(R2>0.5);
for iCell = 1:36
    subplot(6,6,iCell)
    plot(spd_step(2:end), squeeze(spd_resp_all{iarea}(tuned_cells{iarea}(ind(iCell)),2,:))')
    title(num2str(chop(spdSummary(iarea).spdFit(tuned_cells{iarea}(ind(iCell)),2).fit(3:4),2)))
end
    

[p t s] = anovan(C50_sub_all, {area_all});
%% speed preference
fract_topspd = cell(size(area_list));
for iarea = 1:length(area_list)
    expts = unique(expt_id{iarea});
    fract_topspd{iarea} = zeros(1,length(expts));
    [n bin] = histcounts(max_spd{iarea}(tuned_cells{iarea}));
    fprintf([area_list{iarea} '- Pref 30deg/s: ' num2str(n(end)) '/' num2str(sum(n)) '\n'])
    for iexpt = 1:length(expts)
        ind = intersect(tuned_cells{iarea},find(expt_id{iarea} == expts(iexpt)));
        [n bin] = histcounts(max_spd{iarea}(ind));
        fract_topspd{iarea}(1,iexpt) = n(end)./sum(n);
    end
    fprintf([num2str(mean(fract_topspd{iarea},2)) '+/-' num2str(num2str(std(fract_topspd{iarea},[],2)./sqrt(length(fract_topspd{iarea})))) '\n']);
end 

%% fraction responsive
resp_mat = [];
area_mat = [];
spd_mat = [];
spdpref = zeros(length(spd_step),length(area_list));
figure;
for iarea = 1:length(area_list)
    subplot(2,2,1)
    scatter(iarea.*ones(1,length(fract_inc_resp{iarea})), fract_inc_resp{iarea}, 'ob')
    hold on
    errorbar(iarea, mean(fract_inc_resp{iarea},2), std(fract_inc_resp{iarea},[],2)./length(fract_inc_resp{iarea}),'ok')
    ylim([0 1])
    xlim([0 4])
    
    area_list_increspcells{iarea} = [area_list{iarea} '- ' num2str(length(inc_resp_cells{iarea}))];
    
    set(gca, 'XTick', [1:3], 'XTickLabel', area_list)
    subplot(2,2,3)
    resp_mat = [resp_mat; reshape(squeeze(spd_resp_all{iarea}(inc_resp_cells{iarea},2,:)), [length(inc_resp_cells{iarea})*size(spd_resp_all{iarea},3) 1])];
    area_mat = [area_mat; iarea.*ones(length(inc_resp_cells{iarea})*size(spd_resp_all{iarea},3), 1)];
    spd_mat = [spd_mat; reshape(repmat(spd_step(2:end), [length(inc_resp_cells{iarea}) 1]), [length(inc_resp_cells{iarea})*size(spd_resp_all{iarea},3) 1])];
    errorbar(spd_step(2:end),squeeze(mean(spd_resp_all{iarea}(inc_resp_cells{iarea},2,:),1)),  squeeze(std(spd_resp_all{iarea}(inc_resp_cells{iarea},2,:),[],1))./sqrt(length(inc_resp_cells{iarea})))
    hold on
    if iarea==length(area_list)
        legend(area_list_increspcells,'location','northwest')
    end
    subplot(2,2,2)
    scatter(iarea.*ones(1,length(fract_tuned{iarea})), fract_tuned{iarea}, 'ob')
    hold on
    errorbar(iarea, mean(fract_tuned{iarea},2), std(fract_tuned{iarea},[],2)./length(fract_tuned{iarea}),'ok')
    ylim([0 1])
    xlim([0 4])
    area_list_tunedcells{iarea} = [area_list{iarea} '- ' num2str(length(tuned_cells{iarea}))];
    set(gca, 'XTick', [1:3], 'XTickLabel', area_list)
    
    [n bin] = histcounts(max_spd{iarea}(tuned_cells{iarea}));
    spdpref(:,iarea) = n;
end

[p_inc table_inc stats_inc] = anova1(reshape(cell2mat(fract_inc_resp),[length(expt)./length(area_list) length(area_list)]),[],'off');
subplot(2,2,1)
title({'IncrementResp/TotResp',['p = ' num2str(chop(p_inc,2))]})
ylabel('Fraction Inc Responsive')

[p_tune table_tune stats_tune] = anova1(reshape(cell2mat(fract_tuned),[length(expt)./length(area_list) length(area_list)]),[],'off');
subplot(2,2,2)
title({'Tuned/IncrementResp',['p = ' num2str(chop(p_tune,2))]})
ylabel('Fraction Tuned')

subplot(2,2,3)
[p_spdresp table_spdresp stats_spdresp] = anovan(resp_mat, {area_mat,spd_mat});
title('Speed response')
ylabel('dF/F')
xlabel('Speed increment')

subplot(2,2,4)
cats = categorical(cellstr(num2str(spd_step')),cellstr(num2str(spd_step')));
bar(cats,spdpref./sum(spdpref,1))
ylim([0 1])
ylabel('Fraction of cells')
xlabel('Speed step (deg/s)')
legend(area_list_tunedcells,'location','northwest')
title('Speed preferences')
print(fullfile(summary_pn, [mouse_list 'speedSummary_Responsivity.pdf']),'-dpdf','-fillpage')

%kstests for pct correct
[h_PC_V1AL p_PC_V1AL] = kstest2(singleCellPctCorrect{1}(resp_cells{1}),singleCellPctCorrect{2}(resp_cells{2}));
[h_PC_V1PM p_PC_V1PM] = kstest2(singleCellPctCorrect{1}(resp_cells{1}),singleCellPctCorrect{3}(resp_cells{3}));
[h_PC_PMAL p_PC_PMAL] = kstest2(singleCellPctCorrect{3}(resp_cells{3}),singleCellPctCorrect{2}(resp_cells{2}));
%kstests for speed tuning
[h_ST_V1AL p_ST_V1AL] = kstest2(max_spd{1}(tuned_cells{1}),max_spd{2}(tuned_cells{2}));
[h_ST_V1PM p_ST_V1PM] = kstest2(max_spd{1}(tuned_cells{1}),max_spd{3}(tuned_cells{3}));
[h_ST_PMAL p_ST_PMAL] = kstest2(max_spd{3}(tuned_cells{3}),max_spd{2}(tuned_cells{2}));

%% Plots without V1
resp_mat = [];
area_mat = [];
spd_mat = [];
t_area_list = {'AL','PM'};
spdpref = zeros(length(spd_step),length(t_area_list));
area_list_tunedcells = cell(size(t_area_list));
t_area_list_increspcells = cell(size(t_area_list));
figure;
for iarea = 2:3
    t_area_list_increspcells{iarea-1} = [t_area_list{iarea-1} '- ' num2str(length(inc_resp_cells{iarea}))];
    subplot(2,2,1)
    resp_mat = [resp_mat; reshape(squeeze(spd_resp_all{iarea}(inc_resp_cells{iarea},2,:)), [length(inc_resp_cells{iarea})*size(spd_resp_all{iarea},3) 1])];
    area_mat = [area_mat; iarea.*ones(length(inc_resp_cells{iarea})*size(spd_resp_all{iarea},3), 1)];
    spd_mat = [spd_mat; reshape(repmat(spd_step(2:end), [length(inc_resp_cells{iarea}) 1]), [length(inc_resp_cells{iarea})*size(spd_resp_all{iarea},3) 1])];
    errorbar(spd_step(2:end),squeeze(mean(spd_resp_all{iarea}(inc_resp_cells{iarea},2,:),1)),  squeeze(std(spd_resp_all{iarea}(inc_resp_cells{iarea},2,:),[],1))./sqrt(length(inc_resp_cells{iarea})),'-o')
    hold on
    if iarea==3
        legend(t_area_list_increspcells,'location','northwest')
    end
    t_area_list_tunedcells{iarea-1} = [t_area_list{iarea-1} '- ' num2str(length(tuned_cells{iarea}))];
    [n bin] = histcounts(max_spd{iarea}(tuned_cells{iarea}));
    spdpref(:,iarea-1) = n;
end
subplot(2,2,1)
%[p_spdresp table_spdresp stats_spdresp] = anovan(resp_mat, {area_mat,spd_mat});
title('Speed response')
ylabel('dF/F')
xlabel('Speed increment')
ylim([0 0.15])

subplot(2,2,2)
cats = categorical(cellstr(num2str(spd_step')),cellstr(num2str(spd_step')));
bar(cats,spdpref./sum(spdpref,1))
ylim([0 1])
ylabel('Fraction of cells')
xlabel('Speed step (deg/s)')
legend(t_area_list_tunedcells,'location','northwest')
title('Speed preferences')
print(fullfile(summary_pn, [mouse_list 'speedSummary_Responsivity_noV1.pdf']),'-dpdf','-fillpage')

pref_maxspeed=[size(tuned_cells{2},2) spdpref(end,1);size(tuned_cells{3},2) spdpref(end,2)];
 e = sum(pref_maxspeed,2)*sum(pref_maxspeed)/sum(pref_maxspeed(:)); % expected
 X2 = (pref_maxspeed-e).^2./e;
 X2 = sum(X2(:)); % chi square
 df = prod(size(pref_maxspeed)-[1 1]); % degree of freedom
 p_prefmax = 1-chi2cdf(X2,df);
 [AL_m AL_CI] = binofit( spdpref(end,1), size(tuned_cells{2},2));
 [PM_m PM_CI] = binofit( spdpref(end,2), size(tuned_cells{3},2));
 fprintf(['MaxSpeed CI: AL- ' num2str(chop(AL_CI,2)) '; PM- ' num2str(chop(PM_CI,2)) '\n' ]);
fprintf(['MaxSpeed: ' num2str(chop(p_prefmax,2)) '\n']);


R2_AL = [spdSummary(2).spdFit(:,2).Rsq];
R2_PM = [spdSummary(3).spdFit(:,2).Rsq];
R2_fit = [size(tuned_cells{2},2) length(find(R2_AL(tuned_cells{2})>0.5)); size(tuned_cells{3},2) length(find(R2_PM(tuned_cells{3})>0.5))];
 e = sum(R2_fit,2)*sum(R2_fit)/sum(R2_fit(:)); % expected
 X2 = (R2_fit-e).^2./e;
 X2 = sum(X2(:)); % chi square
 df = prod(size(R2_fit)-[1 1]); % degree of freedom
 p_prefmax = 1-chi2cdf(X2,df);
 [AL_m AL_CI] = binofit( length(find(R2_AL(tuned_cells{2})>0.5)), size(tuned_cells{2},2));
 [PM_m PM_CI] = binofit( length(find(R2_PM(tuned_cells{3})>0.5)), size(tuned_cells{3},2));
 fprintf(['R2 fit CI: AL- ' num2str(chop(AL_CI,2)) '; PM- ' num2str(chop(PM_CI,2)) '\n' ]);
fprintf(['R2 fit: ' num2str(chop(p_prefmax,2)) '\n']);

%%
%distribution of percent correct by cell group
figure;
for iarea = 1:length(area_list)
    subplot(2,2,1)
    cdfplot(singleCellPctCorrect{iarea})
    title('All cells')
    area_list_allcells{iarea} = [area_list{iarea} '- ' num2str(length(singleCellPctCorrect{iarea}))];
    if iarea == length(area_list)
        legend(area_list_allcells,'location','northwest')
        vline(0.55)
    end
    hold on
    subplot(2,2,2)
    cdfplot(singleCellPctCorrect{iarea}(resp_cells{iarea}))
    title('Responsive cells')
    area_list_respcells{iarea} = [area_list{iarea} '- ' num2str(length(singleCellPctCorrect{iarea}(resp_cells{iarea})))];
    if iarea == length(area_list)
        legend(area_list_respcells,'location','northwest')
        vline(0.55)
    end
    hold on
    subplot(2,2,3)
    cdfplot(singleCellPctCorrect{iarea}(inc_resp_cells{iarea}))
    title('Increment Responsive cells')
    area_list_increspcells{iarea} = [area_list{iarea} '- ' num2str(length(singleCellPctCorrect{iarea}(inc_resp_cells{iarea})))];
    if iarea == length(area_list)
        legend(area_list_increspcells,'location','northwest')
        vline(0.55)
    end
    hold on
    subplot(2,2,4)
    cdfplot(singleCellPctCorrect{iarea}(tuned_cells{iarea}))
    title('Tuned cells')
    area_list_tunedcells{iarea} = [area_list{iarea} '- ' num2str(length(singleCellPctCorrect{iarea}(tuned_cells{iarea})))];
    if iarea == length(area_list)
        legend(area_list_tunedcells,'location','northwest')
        vline(0.55)
    end
    hold on
end

print(fullfile(summary_pn, [mouse_list 'speedSummary_singleCellDecodeCDF.pdf']),'-dpdf','-bestfit')


figure;
for iarea = 1:length(area_list)
    subplot(2,2,iarea)
    scatter(spd_step(max_spd{iarea}+1), singleCellPctCorrect{iarea});
    ylim([0 1])
    xlim([0 35])
    ylabel('% correct decoding')
    xlabel('Speed step (deg/s)')
    title(area_list{iarea})
    hold on
    scatter(spd_step(max_spd{iarea}(resp_cells{iarea})+1), singleCellPctCorrect{iarea}(resp_cells{iarea}));
    scatter(spd_step(max_spd{iarea}(tuned_cells{iarea})+1), singleCellPctCorrect{iarea}(tuned_cells{iarea}));
    if iarea == 1
        legend({'all', 'responsive', 'tuned'},'location','southeast')
    end
end

print(fullfile(summary_pn, [mouse_list 'speedSummary_singleCellDecodeBySpeedPref.pdf']),'-dpdf','-bestfit')

figure;
for iarea = 1:length(area_list)
    subplot(2,2,iarea)
    scatter(spd_step(max_spd{iarea}+1), cell_weights{iarea});
    xlim([0 35])
    ylim([-20 100])
    ylabel('Weight')
    xlabel('Speed step (deg/s)')
    title(area_list{iarea})
    hold on
    scatter(spd_step(max_spd{iarea}(resp_cells{iarea})+1), cell_weights{iarea}(resp_cells{iarea}));
    scatter(spd_step(max_spd{iarea}(tuned_cells{iarea})+1), cell_weights{iarea}(tuned_cells{iarea}));
    if iarea == 3
        legend({'all', 'responsive', 'tuned'},'location','northwest')
    end
    subplot(2,2,4)
    scatter(singleCellPctCorrect{iarea}, cell_weights{iarea});
    xlim([0 1])
    ylim([-20 100])
    ylabel('Weight')
    xlabel('% correct decoding')
    hold on
    if iarea == 3
        legend(area_list,'location','northwest')
    end
end

print(fullfile(summary_pn, [mouse_list 'speedSummary_singleCellWeightBySpeedPref.pdf']),'-dpdf','-bestfit')

figure;
subplot(2,3,1)
for iarea = 1:length(area_list)
    cdfplot(spd_step(max_spd{iarea}(inc_resp_cells{iarea})+1));
    ylim([0 1])
    xlim([0 35])
    ylabel('F(x)')
    xlabel('Speed step (deg/s)')
    hold on
end
title('Dots- steps from base')
subplot(2,3,2)
for iarea = 1:length(area_list)
    cdfplot(spd_step(max_spd_nobase{iarea}(resp_cells{iarea})+1));
    ylim([0 1])
    xlim([0 35])
    ylabel('F(x)')
    xlabel('Speed step (deg/s)')
    hold on
end
title('Dots- no base')
subplot(2,3,3)
for iarea = 1:length(area_list)
    cdfplot(grating_spds(max_spd_grating{iarea}(TF_resp_cells{iarea})));
    ylim([0 1])
    xlim([0 70])
    ylabel('F(x)')
    xlabel('Speed step (deg/s)')
    hold on
end
title('Gratings- no base')
subplot(2,3,4)
for iarea = 1:length(area_list)
    cdfplot(spd_step(max_spd{iarea}(tuned_cells{iarea})+1));
    ylim([0 1])
    xlim([0 35])
    ylabel('F(x)')
    xlabel('Speed step (deg/s)')
    hold on
end
subplot(2,3,5)
for iarea = 1:length(area_list)
    cdfplot(spd_step(max_spd_nobase{iarea}(tuned_cells{iarea})+1));
    ylim([0 1])
    xlim([0 35])
    ylabel('F(x)')
    xlabel('Speed step (deg/s)')
    hold on
end
subplot(2,3,6)
for iarea = 1:length(area_list)
    cdfplot(grating_spds(max_spd_grating{iarea}(TF_tuned_cells{iarea})));
    ylim([0 1])
    xlim([0 70])
    ylabel('F(x)')
    xlabel('Speed step (deg/s)')
    hold on
end
legend(area_list, 'location','southeast')
print(fullfile(summary_pn, [mouse_list 'speedSummary_speedPrefCDF.pdf']),'-dpdf','-bestfit')

