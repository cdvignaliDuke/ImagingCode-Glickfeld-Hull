clear all;
doRedChannel = 0;
ds = 'CrossOriRandDir_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 15;
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandDirSummary');

stim_OSI_all = [];
plaid_OSI_all = [];
stim_DSI_all = [];
plaid_DSI_all = [];
Rc_all = [];
Rp_all = [];
plaid_SI_all = [];
totCells = 0;
resp_ind_all = [];
resp_ind_dir_all = [];
resp_ind_plaid_all = [];

for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    fprintf([mouse ' ' date '\n'])

    %% load data

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dirAnalysis.mat']), 'Rp', 'Rc', 'stim_OSI', 'stim_DSI', 'plaid_OSI', 'plaid_DSI', 'plaid_SI', 'nCells');

    stim_OSI_all = [stim_OSI_all stim_OSI];
    plaid_OSI_all = [plaid_OSI_all plaid_OSI];
    stim_DSI_all = [stim_DSI_all stim_DSI];
    plaid_DSI_all = [plaid_DSI_all plaid_DSI];
    Rc_all = [Rc_all Rc];
    Rp_all = [Rp_all Rp];
    plaid_SI_all = [plaid_SI_all plaid_SI];
    
    resp_ind = find(sum(sum(h_resp,2),3));
    resp_ind_dir = find(sum(h_resp(:,:,1),2));
    resp_ind_plaid = find(sum(h_resp(:,:,2),2));
    
    resp_ind_all = [resp_ind_all resp_ind'+totCells];
    resp_ind_dir_all = [resp_ind_dir_all resp_ind_dir'+totCells];
    resp_ind_plaid_all = [resp_ind_plaid_all resp_ind_plaid'+totCells];
    
    totCells = totCells+nCells;
    

end


figure;
subplot(2,2,1)
cdfplot(stim_OSI_all(resp_ind_all))
hold on
cdfplot(plaid_OSI_all(resp_ind_all))
xlabel('OSI')
legend({'Stim','Plaid'},'Location','southeast')
title('')
subplot(2,2,2)
cdfplot(stim_DSI_all(resp_ind_all))
hold on
cdfplot(plaid_DSI_all(resp_ind_all))
xlabel('DSI')
title('')
legend({'Stim','Plaid'},'Location','southeast')
subplot(2,2,3)
cdfplot(Rc_all(resp_ind_all))
hold on
cdfplot(Rp_all(resp_ind_all))
xlabel('Rc/Rp')
xlim([-2 10])
title('')
legend({'Rc','Rp'},'Location','southeast')
subplot(2,2,4)
cdfplot(plaid_SI_all(resp_ind_all))
hold on
cdfplot(plaid_SI_all(intersect(resp_ind_all,find(stim_OSI_all<0.5))))
cdfplot(plaid_SI_all(intersect(resp_ind_all,find(stim_DSI_all<0.5))))
xlabel('Suppression Index')
title('')
legend({'All','stim OSI<0.5', 'stim DSI<0.5'},'Location','southeast')
suptitle(['All responsive cells- n = ' num2str(length(resp_ind_all))])
print(fullfile(summaryDir, 'randDir_OSI-DSI-Rc-Rp-SI_Summary.pdf'),'-dpdf', '-fillpage')       


figure;
scatter(Rc_all(resp_ind_plaid_all),Rp_all(resp_ind_plaid_all))
xlabel('Rc')
ylabel('Rp')
xlim([-2 20])
ylim([-2 20])
axis square
title(['Plaid responsive cells- n = ' num2str(length(resp_ind_plaid_all))])
print(fullfile(summaryDir, 'randDir_Rc-Rp_Scatter.pdf'),'-dpdf', '-fillpage') 

figure;
subplot(3,2,1)
cdfplot(stim_OSI_all(intersect(resp_ind_all,find(plaid_SI_all<0))))
hold on
cdfplot(stim_OSI_all(intersect(resp_ind_all,find(plaid_SI_all>0))))
xlabel('stim OSI')
legend({'SI<0', 'SI>0'},'Location','northwest')
title('')
subplot(3,2,2)
cdfplot(plaid_OSI_all(intersect(resp_ind_all,find(plaid_SI_all<0))))
hold on
cdfplot(plaid_OSI_all(intersect(resp_ind_all,find(plaid_SI_all>0))))
xlabel('plaid OSI')
title('')
subplot(3,2,3)
cdfplot(stim_DSI_all(intersect(resp_ind_all,find(plaid_SI_all<0))))
hold on
cdfplot(stim_DSI_all(intersect(resp_ind_all,find(plaid_SI_all>0))))
xlabel('stim DSI')
title('')
subplot(3,2,4)
cdfplot(plaid_DSI_all(intersect(resp_ind_all,find(plaid_SI_all<0))))
hold on
cdfplot(plaid_DSI_all(intersect(resp_ind_all,find(plaid_SI_all>0))))
xlabel('plaid DSI')
title('')
subplot(3,2,5)
cdfplot(Rc_all(intersect(resp_ind_all,find(plaid_SI_all<0))))
hold on
cdfplot(Rc_all(intersect(resp_ind_all,find(plaid_SI_all>0))))
xlabel('Rc')
xlim([-2 10])
title('')
subplot(3,2,6)
cdfplot(Rp_all(intersect(resp_ind_all,find(plaid_SI_all<0))))
hold on
cdfplot(Rp_all(intersect(resp_ind_all,find(plaid_SI_all>0))))
xlabel('Rp')
xlim([-2 10])
title('')
suptitle({'High vs low Suppression index',['All responsive cells- n = ' num2str(length(resp_ind_all))]})
print(fullfile(summaryDir, 'randDir_highVlowSI.pdf'),'-dpdf', '-fillpage') 

figure;
subplot(2,2,1)
cdfplot(plaid_SI_all(intersect(resp_ind_all,find(stim_OSI_all<0.5))))
hold on
cdfplot(plaid_SI_all(intersect(resp_ind_all,find(stim_OSI_all>0.5))))
xlabel('Suppression index')
legend({'OSI<0.5', 'OSI>0.5'},'Location','northwest')
title('')
subplot(2,2,2)
cdfplot(stim_DSI_all(intersect(resp_ind_all,find(stim_OSI_all<0.5))))
hold on
cdfplot(stim_DSI_all(intersect(resp_ind_all,find(stim_OSI_all>0.5))))
xlabel('stim DSI')
title('')
subplot(2,2,3)
cdfplot(Rc_all(intersect(resp_ind_all,find(stim_OSI_all<0.5))))
hold on
cdfplot(Rc_all(intersect(resp_ind_all,find(stim_OSI_all>0.5))))
xlabel('Rc')
xlim([-2 10])
title('')
subplot(2,2,4)
cdfplot(Rp_all(intersect(resp_ind_all,find(stim_OSI_all<0.5))))
hold on
cdfplot(Rp_all(intersect(resp_ind_all,find(stim_OSI_all>0.5))))
xlabel('Rp')
xlim([-2 10])
title('')
suptitle({'High vs low OSI', ['All responsive cells- n = ' num2str(length(resp_ind_all))]})
print(fullfile(summaryDir, 'randDir_highVlowOSI.pdf'),'-dpdf', '-fillpage') 

figure;
subplot(2,2,1)
cdfplot(plaid_SI_all(intersect(resp_ind_all,find(stim_DSI_all<0.5))))
hold on
cdfplot(plaid_SI_all(intersect(resp_ind_all,find(stim_DSI_all>0.5))))
xlabel('Suppression index')
legend({'DSI<0.5', 'DSI>0.5'},'Location','northwest')
title('')
subplot(2,2,2)
cdfplot(plaid_DSI_all(intersect(resp_ind_all,find(stim_DSI_all<0.5))))
hold on
cdfplot(plaid_DSI_all(intersect(resp_ind_all,find(stim_DSI_all>0.5))))
xlabel('plaid DSI')
title('')
subplot(2,2,3)
cdfplot(Rc_all(intersect(resp_ind_all,find(stim_DSI_all<0.5))))
hold on
cdfplot(Rc_all(intersect(resp_ind_all,find(stim_DSI_all>0.5))))
xlabel('Rc')
xlim([-2 10])
title('')
subplot(2,2,4)
cdfplot(Rp_all(intersect(resp_ind_all,find(stim_DSI_all<0.5))))
hold on
cdfplot(Rp_all(intersect(resp_ind_all,find(stim_DSI_all>0.5))))
xlabel('Rp')
xlim([-2 10])
title('')
suptitle({'High vs low DSI', ['All responsive cells- n = ' num2str(length(resp_ind_all))]})
print(fullfile(summaryDir, 'randDir_highVlowDSI.pdf'),'-dpdf', '-fillpage') 
