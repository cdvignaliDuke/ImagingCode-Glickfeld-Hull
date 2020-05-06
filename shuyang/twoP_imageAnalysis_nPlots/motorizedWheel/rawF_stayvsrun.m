% compare the raw fluorescence during running vs. stationary WHEN THE CELL
% IS NOT FIRING. 
% doing this analysis because the baseline of CaAmp plot is different for
% stationary and running and I need to figure out why

clear;
sessions = {'200116_img1041','200214_img1042','200217_img1061','200225_img1049','200319_img1064','200319_img1064_2'};
days = {'1041-200116_1','1042-200214_1','1061-200217_1','1049-200225_1','1064-200319_1','1064-200319_2'};
%sessions = {'200116_img1041','200225_img1049','200319_img1064_2'};
%days = {'1041-200116_1','1049-200225_1','1064-200319_2'};
image_analysis_base = 'Z:\Analysis\motorizedWheel_Analysis\running\imaging_analysis\'; 
color_code = {'c','r','y','g'};

%% 
for ii = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    rawF_output = load([image_analysis_dest sessions{ii},'_deconvolution_thresh-4_TCave_cl.mat']);
    rawF = rawF_output.TCave_cl;
    CaAmp_output = load([image_analysis_dest sessions{ii} '_isoCaEvent.mat']);
    stay_CSpinx = CaAmp_output.spk_iso_stay_neurons; % when peak of Ca event happens
    runfast_CSpinx = CaAmp_output.spk_iso_runfast_neurons;
    runslow_CSpinx = CaAmp_output.spk_iso_runslow_neurons;
    
    stay_befoCS = cell(size(stay_CSpinx));
    runfast_befoCS = cell(size(stay_CSpinx));
    runslow_befoCS = cell(size(stay_CSpinx));
    rawF_stay_befoCS = cell(size(stay_CSpinx));
    rawF_runfast_befoCS = cell(size(stay_CSpinx));
    rawF_runslow_befoCS = cell(size(stay_CSpinx));
    rawF_stay_befoCS_vec = zeros(1,size(stay_CSpinx,2));
    rawF_runfast_befoCS_vec = zeros(1,size(stay_CSpinx,2));
    rawF_runslow_befoCS_vec = zeros(1,size(stay_CSpinx,2));
    
    for c = 1:size(stay_CSpinx,2)
        for s1 = 1:length(stay_CSpinx{c}) %for each spike, take the baseline period before Ca transient
            stay_befoCS{c} = [stay_befoCS{c} stay_CSpinx{c}(s1)-18:stay_CSpinx{c}(s1)-4]; %peaks are at least 18 frames from each other (600ms, each frame is 33Hz.), and there're about 4 raising frames
        end
        
        for s2 = 1:length(runfast_CSpinx{c}) %for each spike, take the baseline period before Ca transient
            runfast_befoCS{c} = [runfast_befoCS{c} runfast_CSpinx{c}(s2)-18:runfast_CSpinx{c}(s2)-4];
        end
        for s3 = 1:length(runslow_CSpinx{c}) %for each spike, take the baseline period before Ca transient
            runslow_befoCS{c} = [runslow_befoCS{c} runslow_CSpinx{c}(s3)-18:runslow_CSpinx{c}(s3)-4];
        end
        
        
        rawF_stay_befoCS{c} = rawF(stay_befoCS{c},c);
        rawF_runfast_befoCS{c} = rawF(runfast_befoCS{c},c);
        rawF_runslow_befoCS{c} = rawF(runslow_befoCS{c},c);
        
        rawF_stay_befoCS_vec(c) = mean(rawF_stay_befoCS{c}); %for each cell, average across time
        rawF_runfast_befoCS_vec(c) = mean(rawF_runfast_befoCS{c}); %for each cell, average across time
        rawF_runslow_befoCS_vec(c) = mean(rawF_runslow_befoCS{c}); %for each cell, average across time    
    end
    
    rawF_compare = [rawF_stay_befoCS_vec',rawF_runfast_befoCS_vec',rawF_runslow_befoCS_vec'];
    rawF_compare(~all(~isnan(rawF_compare),2),:)=[];
    %not all of the cells fire during slow running and fast running, so there will be lines in rawF_compare has NaN. delete those lines/cells
    
    
    x = [1,2,3];
    x_plot = repmat(x,size(rawF_compare,1),1);
    figure; 
    plot(x_plot',rawF_compare','.','LineStyle','-','linewidth', 1.25,'MarkerSize',20,'color',[0.5294 0.5294 0.5294]); hold on;
    fast_errbar(x,rawF_compare,1,'color',[0.8896 0.1922 0.1557]);
    xlabel ('behavioral state');xlim([0.5 3.5]); axis square;
    set(gca,'XTick',x,'XTicklabel',{'stationary','fast running','slow running'});
    ylabel('raw F w/o CS');
    title(sessions{ii});
    savefig([image_analysis_dest sessions{ii} '_rawF_baseline.fig']);
    save([image_analysis_dest sessions{ii},'_deconvolution_thresh-4_TCave_cl.mat'],'rawF_compare','-append');
end


%% across sessions
rawF_compare_all = [];

for ii = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    rawF_output = load([image_analysis_dest sessions{ii},'_deconvolution_thresh-4_TCave_cl.mat']);
    rawF_compare = rawF_output.rawF_compare;
    rawF_compare_all = cat(1,rawF_compare_all,rawF_compare);
end

x = [1,2,3];
figure;
x_plot = repmat(x,size(rawF_compare_all,1),1);
plot(x_plot',rawF_compare_all','.','LineStyle','-','linewidth', 1.25,'MarkerSize',20,'color',[0.5294 0.5294 0.5294]); hold on;
fast_errbar(x,rawF_compare_all,1,'color',[0.8896 0.1922 0.1557]);
xlabel ('behavioral state');xlim([0.5 3.5]); axis square;
set(gca,'XTick',x,'XTicklabel',{'stationary','fast running','slow running'});
ylabel('raw F w/o CS');
title('across sessions');
savefig(['Z:\Analysis\motorizedWheel_Analysis\running\across_sessions' '_rawF_baseline.fig']);



