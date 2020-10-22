clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_lowSF_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 30;
nexp = size(expt,2);

for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

    fprintf(['2P imaging sine fitting analysis\nSelected data:\nMouse: ' mouse '\nDate: ' date '\nExperiments:\n'])
    for irun=1:nrun
        fprintf([ImgFolder{irun} ' - ' time{irun} '\n'])
    end


    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))
    
    rAM = data_dfof_con_ph_tc_avg(:,resp_ind,1,2,1);
    nframes = size(rAM,1);
    rBM = cell(1,nMaskPhas);
    figure; 
    for i = 1:nMaskPhas
        rBM{i} = data_dfof_con_ph_tc_avg(:,resp_ind,3,2,i);
        subplot(2,2,i)
        scatter(mean(rAM,1),mean(rBM{i},1));
        title(['Phase = ' num2str(maskPhas(i))])
        xlim([0 1.5])
        xlabel('Stim')
        ylim([0 1.5])
        xlabel('Stim + Mask')
    end
    
    figure; 
    plot(mean(rAM,2))
    X = smoothdata(rAM,1);
    legstr = 'grating';
    hold on
    for i = 1:nMaskPhas
        plot(mean(rBM{i},2))
        legstr = strvcat(legstr, [num2str(maskPhas(i)) ' Mask']);
        X = [X; smoothdata(rBM{i})];
    end
    legend(legstr, 'location', 'northwest')
    
    figure;
    [Y,e]=cmdscale(pdist(X,'cos'));
    plot(e(1:20),'k.-','markersize',25)
    box off
    xlabel('Rank')
    ylabel('Eigenvalue')
    
    figure;
    u = 1:nframes;   % indices for masked and unmasked conditions
    yA = Y(u,:);
    yB = cell(1,nMaskPhas);
    m = u;
    for i = 1:nMaskPhas
        subplot(2,2,i)
        plot3(yA(:,1),yA(:,2),yA(:,3),'b-','LineWidth',2)
        hold on
        m = m+nframes;
        yB{i} = Y(m,:);
        plot3(yB{i}(:,1),yB{i}(:,2),yB{i}(:,3),'r-','LineWidth',2)
        title(num2str(maskPhas(i)))
    end
        view([7.244 11.396])
        hold off
    
    figure;
    plot3(yA(:,1),yA(:,2),yA(:,3),'LineWidth',2)
    legstr = cell(1,nMaskPhas+1);
    legstr{1} = 'Grating';
    hold on
    for i = 1:nMaskPhas
        plot3(yB{i}(:,1),yB{i}(:,2),yB{i}(:,3),'LineWidth',2)
        hold on
        legstr{i+1} = num2str(maskPhas(i));
    end
    view([7.244 11.396])
    hold off
    legend(legstr)
        
    ndim = 6;            % select embedding dimension
    figure
    for i = 1:nMaskPhas
        yA = yA(:,1:ndim);
        yB{i} = yB{i}(:,1:ndim);

        yAn = bsxfun(@rdivide,yA,sqrt(sum(yA.^2,2))); % normalize responses
        yBn = bsxfun(@rdivide,yB{i},sqrt(sum(yB{i}.^2,2)));

        yAh = [yAn ones(size(yA,1),1)]; % put in homogeneous coordinates 
        yBh = [yBn ones(size(yB{i},1),1)];

        A = pinv(yAh)*yBh;              % estimate affine transform
        yBhp = yAh*A;                   % estimate prediction
        subplot(2,2,i)
        plot3(yAh(:,1),yAh(:,2),yAh(:,3),'b-','LineWidth',2)
        hold on
        plot3(yBh(:,1),yBh(:,2),yBh(:,3),'r-','LineWidth',2)
        plot3(yBhp(:,1),yBhp(:,2),yBhp(:,3),'g-','LineWidth',2)
        axis equal
        title(num2str(maskPhas(i)))
        view([-149.29 12.17])
    end
    
end 
