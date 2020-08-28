clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRotate_ExptList';
rc = behavConstsAV;
eval(ds)
nexp = size(expt,2);
frame_rate = 15;

yA_all = cell(2,1);
yB_all = cell(2,4);
for ir = 1:2
    yA_all{ir} = [];
    for ip = 1:4
        yB_all{ir,ip} = [];
    end
end
nCells_all = zeros(nexp,1);
fprintf('Cross ori rotate summary \n')

for iexp = 1:nexp

    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    ImgFolder = expt(iexp).coFolder;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

    fprintf([mouse ' ' date '\n'])

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dataStim.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_TCs.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_input.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respAvg.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_cosSpace.mat']))
    
    nCells = size(npSub_tc,2);
    nCells_all(iexp,1) = nCells;
    if length(yA)<2
        yA{2} = nan(size(yA{1}));
        for ip = 1:nMaskPhas
            yB{2,ip} = nan(size(yB{1,ip}));
        end
    end
    for ir = 1:2
        yA_all{ir} = cat(3, yA_all{ir}, yA{ir}(:,1:100));
        for ip = 1:nMaskPhas
            yB_all{ir,ip} = cat(3, yB_all{ir,ip}, yB{ir,ip}(:,1:100));
        end
    end
end
%%
Summary_fn = fullfile(LG_base,  'Analysis\2P', 'CrossOri', 'RotateSummary');
mask_str = {'stationary','rotate'};
yA_all_avg = cell(size(yA_all));
yB_all_avg = cell(size(yB_all));
for ir = 1:nR
    figure;
    for iexp = 1:nexp
        subplot(3,2,iexp)
        plot(yA_all{ir}(:,1,iexp),yA_all{ir}(:,2,iexp),'LineWidth',2)
        hold on
        legstr = cell(1,nMaskPhas+1);
        legstr{1} = 'Grating';
        for i = 1:nMaskPhas
            plot(yB_all{ir,i}(:,1,iexp),yB_all{ir,i}(:,2,iexp),'LineWidth',2)
            hold on
            legstr{i+1} = num2str(maskPhas(i));
        end
        
        title([expt(iexp).mouse ' ' expt(iexp).date '- ' num2str(nCells_all(iexp,:)) ' cells'])
        legend(legstr)
    end
    
    movegui('center')
    suptitle(['Mask- ' mask_str(ir)])
    yA_all_avg{ir} = nanmean(yA_all{ir},3);
    for i = 1:nMaskPhas
        yB_all_avg{ir,i} = nanmean(yB_all{ir,i},3);
    end
    print(fullfile(Summary_fn, ['allexpt_mask' mask_str{ir} '_2D.pdf']),'-dpdf','-fillpage')
end

for ir = 1:nR
    figure;
    for iexp = 1:nexp
        subplot(3,2,iexp)
        plot3(yA_all{ir}(:,1,iexp),yA_all{ir}(:,2,iexp),yA_all{ir}(:,3,iexp),'LineWidth',2)
        hold on
        legstr = cell(1,nMaskPhas+1);
        legstr{1} = 'Grating';
        for i = 1:nMaskPhas
            plot3(yB_all{ir,i}(:,1,iexp),yB_all{ir,i}(:,2,iexp),yB_all{ir,i}(:,3,iexp),'LineWidth',2)
            hold on
            legstr{i+1} = num2str(maskPhas(i));
        end
        
        title([expt(iexp).mouse ' ' expt(iexp).date '- ' num2str(nCells_all(iexp,:)) ' cells'])
        legend(legstr)
        view([7.244 11.396])
    end
    
    movegui('center')
    suptitle(['Mask- ' mask_str(ir)])
    print(fullfile(Summary_fn, ['allexpt_mask' mask_str{ir} '_3D.pdf']),'-dpdf','-fillpage')
end

for ir = 1:nR
    figure;
    plot(yA_all_avg{ir}(:,1),yA_all_avg{ir}(:,2),'LineWidth',2)
    hold on
    legstr = cell(1,nMaskPhas+1);
    legstr{1} = 'Grating';
    for i = 1:nMaskPhas
         plot(yB_all_avg{ir,i}(:,1),yB_all_avg{ir,i}(:,2),'LineWidth',2)
        hold on
        legstr{i+1} = num2str(maskPhas(i));
    end
    movegui('center')
title(['Mask- ' mask_str(ir)])
legend(legstr)
print(fullfile(Summary_fn, ['allexpt_mask' mask_str{ir} '_2D_avg.pdf']),'-dpdf','-fillpage')
end
