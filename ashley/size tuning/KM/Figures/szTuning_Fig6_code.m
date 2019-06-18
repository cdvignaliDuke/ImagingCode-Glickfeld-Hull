rc = behavConstsAV;
mouse = '1205';
expDate = '190131';
dataFolder = '003';
i = 164;

fnout = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging',expDate);
load(fullfile(fnout, dataFolder, [mouse '_' expDate '_input.mat']))
load(fullfile(fnout, dataFolder, [mouse '_' expDate '_lbub_fits.mat']))
load(fullfile(fnout, dataFolder, [mouse '_' expDate '_Fit_struct_sub.mat']))
load(fullfile(fnout, dataFolder, [mouse '_' expDate '_Tuning.mat']))
%%
sz_mat = celleqel2mat_padded(input.tGratingDiameterDeg);
szs = unique(sz_mat);
szRng = linspace(0,max(szs));
sizeMean = squeeze(tuning_mat(:,1,:));
sizeSEM = squeeze(tuning_mat(:,2,:));
Fit_struct = Fit_struct_sub;
clear Fit_struct_sub

figure;
s = Fit_struct(i).True.s_;
errorbar([0 szs],[0 sizeMean(:,i)'],[0 sizeSEM(:,i)'],'ok')
hold on
plot(s.szs0,s.data,'.b')
plot(szRng,s.fitout1,'-r')
plot(szRng,s.fitout2,'-g')
hold off
ylim([min([-0.5*s.maxResp1 min(s.data)]) 1.2*max([s.maxResp2 max(s.data)])])
title([mouse ' ' expDate '- Cell #' num2str(i)]);
xlabel('Stimulus size (deg)')
ylabel('dF/F')
ylim([-0.15 0.4])
axis square
print(['\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Manuscripts\sizeTuning\Figures5-6\i1205_190131_Cell#' num2str(i) '_sizeTuning.pdf'],'-dpdf','-bestfit')
