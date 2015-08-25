%% find which cells resp to which flash
baselineITI_V = [];
prevflashbaseline_V = [];
flashresp_V = [];
baselineITI_AV = [];
prevflashbaseline_AV = [];
flashresp_AV = [];

flashresp2ITI_V = zeros(nCells,length(cycles));
flashresp2prevflash_V = zeros(nCells,length(cycles));
flashresp2ITI_AV = zeros(nCells,length(cycles));
flashresp2prevflash_AV = zeros(nCells,length(cycles));


for icyc = 1:length(cycles)
    tempdata = cycDataDFoverF_cmlvNoTarget{icyc};
    v_ind = cycV_ind{icyc};
    av_ind = cycAV_ind{icyc};
    baselineITI_V{icyc} = squeeze(mean(tempdata(1:30,:,v_ind),1));
    baselineITI_AV{icyc} = squeeze(mean(tempdata(1:30,:,av_ind),1));
    prevflashbaseline_V{icyc} = squeeze(mean(tempdata( (28+(cycTime*(cycles(icyc)-1))):(28+(cycTime*(cycles(icyc)-1)))+4 ,:,v_ind),1));
    prevflashbaseline_AV{icyc} = squeeze(mean(tempdata( (28+(cycTime*(cycles(icyc)-1))):(28+(cycTime*(cycles(icyc)-1)))+4 ,:,av_ind),1));
    flashresp_V{icyc} = squeeze(mean(tempdata( (33+(cycTime*(cycles(icyc)-1))):(33+(cycTime*(cycles(icyc)-1)))+4 ,:,v_ind),1));
    flashresp_AV{icyc} = squeeze(mean(tempdata( (33+(cycTime*(cycles(icyc)-1))):(33+(cycTime*(cycles(icyc)-1)))+4 ,:,av_ind),1));
    
    flashresp_V_mean(:,icyc) = mean(flashresp_V{icyc},2);
    flashresp_V_var(:,icyc) = var(flashresp_V{icyc},[],2);
    flashresp_AV_mean(:,icyc) = mean(flashresp_AV{icyc},2);
    flashresp_AV_var(:,icyc) = var(flashresp_AV{icyc},[],2);
    
    flashresp2ITI_V(:,icyc) = ttest(baselineITI_V{icyc}',flashresp_V{icyc}','alpha', 0.01);
    flashresp2prevflash_V(:,icyc) = ttest(prevflashbaseline_V{icyc}',flashresp_V{icyc}','alpha', 0.01);
    flashresp2ITI_AV(:,icyc) = ttest(baselineITI_AV{icyc}',flashresp_AV{icyc}','alpha', 0.01);
    flashresp2prevflash_AV(:,icyc) = ttest(prevflashbaseline_AV{icyc}',flashresp_AV{icyc}','alpha', 0.01);
    
    flashrespthresh_V{icyc} = find(any(mean(tempdata(:,:,v_ind),3) >0.05,1));
    flashrespthresh_AV{icyc} = find(any(mean(tempdata(:,:,av_ind),3) >0.05,1));
end

baseflashrespcells_V = find(any(flashresp2ITI_V == 1,2));
baseflashrespcells_AV = find(any(flashresp2ITI_AV == 1,2));

%%
figure;
scatter(flashresp_V_mean(:,7),flashresp_V_var(:,7),'g')
hold on
scatter(flashresp_AV_mean(:,7),flashresp_AV_var(:,7),'k')
hold on
refline(1,0)
hold on
% xlim([-0.08 0.14]);
% ylim([-0.08 0.14]);

%% plot driven cells
CYC = 2;
drivencells = flashrespthresh_V{CYC};
tempdata = cycDataDFoverF_cmlvNoTarget{CYC};
tempVind = cycV_ind{CYC};
tempAVind = cycAV_ind{CYC};

for i = 1:length(drivencells)
figure;
% 
%     subplot(5,5,i)
plot(tsmovavg(squeeze(tempdata(:,(drivencells(i)),tempVind)),'s',3,1),'g')
alpha(0.25)
hold on
plot(tsmovavg(squeeze(tempdata(:,(drivencells(i)),tempAVind)),'s',3,1),'k')
% hold on
vline(30,'k') 
    for i = 1:cycles(CYC)-1
        L = (CYC*cycTime)+31;
        vline(L,'k:');
        hold on
    end
end



%% color-code cells by orientation preference
DirFolder = '006';
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, DirFolder);
cd(fileSave);
% load('oriTuningPreferences.mat')
load('TuningPreferences.mat')

dataTrialStart = cycDataDFoverF_cmlvNoTarget{4};
v_ind = cycV_ind{4};
% a_ind = cycA_ind{1};
a_ind = cycAV_ind{4};


preStimResp_V = zeros(size(v_ind,2),size(dataTrialStart,2));
for itrial =1:size(v_ind,1);
    for icell = 1:size(dataTrialStart,2)
        preStimResp_V(itrial,icell) = mean(dataTrialStart(1:30,icell,v_ind(itrial)),1);
    end
end

baselineStimResp_V = zeros(size(v_ind,2),size(dataTrialStart,2));
for itrial = 1:size(v_ind,1);
    for icell = 1:size(dataTrialStart,2)
        baselineStimResp_V(itrial,icell) = mean(dataTrialStart(36:end,icell,v_ind(itrial)),1);
    end
end

baselineStimRespTtest_V= ttest(preStimResp_V,baselineStimResp_V,'alpha', 0.01);
baselineStimRespIndex_V = find(baselineStimRespTtest_V == 1);


cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
cellsSelectZero_dir = intersect(dirSlctvCells,cellsPrefZero);
cellsSelectZero_ori = intersect(oriSlctvCells,cellsPrefZero);
cellsSelectZero = union(cellsSelectZero_dir,cellsSelectZero_ori);
cellsPrefRespZero = intersect(baselineStimRespIndex_V,cellsPrefZero);

cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);
cellsSelectNinety_dir = intersect(dirSlctvCells,cellsPrefNinety);
cellsSelectNinety_ori = intersect(oriSlctvCells,cellsPrefNinety);
cellsSelectNinety = union(cellsSelectNinety_dir,cellsSelectNinety_ori);
cellsPrefRespNinety = intersect(baselineStimRespIndex_V,cellsPrefNinety);

nCells = size(cycDataDFoverF_cmlvNoTarget{7},2);
oriSlctvCellsAll = union(oriSlctvCells,dirSlctvCells);
notSlctvCells = setdiff([1:nCells],oriSlctvCellsAll);
notRespCells = setdiff([1:nCells],baselineStimRespIndex_V);