%% maps
% El Az, sigEl sigAz, RFsigma (geomean sigEl+sigAz), RF-stim, prefSize (0.8), SI (0.8)

% need to load each exp: masks, ret data, sizeFits
% only use goodfit_size cells

%% load csv to define exp

clear all; clc;
%expfile = '\\CRASH.dhe.duke.edu\data\home\kevin\Code\Ai9x_experiment_list.txt';
%expfile = 'C:\Users\kevin\Documents\Repositories\ImagingCode-Glickfeld-Hull\kevin\Size Tuning Code\Ai9x_experiment_list.txt';
expfile = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Code\Ai9x_experiment_list.txt';
fID = fopen(expfile);
head = textscan(fID,'%s%s%s%s%s',1,'delimiter',',');
head = vertcat(head{:});
temp = textscan(fID,'%s%s%s%s%s','delimiter',',','HeaderLines',1);
temp = horzcat(temp{:});
expdata = cell2table(temp,'VariableNames',head);
nExp = size(expdata,1);
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);

fprintf(['Cell maps, prefSize v Az analysis - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

%% load each exp individually

close all
for iExp = 49%1:nExp
    date = expdata.date{iExp};
    mouse = expdata.mouse{iExp};
    ret_str = 'runs-002';
    run_str = expdata.run_str{iExp};
    fprintf('Experiment %d/%d: mouse %s date %s\n',iExp,nExp,mouse,date)
    
    % masks
    load(fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_mask_cell.mat']))
    fprintf('Loaded: cell masks')
    
    % lbub_fits (ret)
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' ret_str], [date '_' mouse '_' ret_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' ret_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind')
    nCellsR = size(lbub_fits,1);
    lbub_fits_ret = lbub_fits;
    
    %cellAz_norm = lbub_fits(goodfit_ind,4,4) - stimPosExp(i,2);
    %cellEl_norm = lbub_fits(goodfit_ind,5,4) - stimPosExp(i,1);
    sigmax = lbub_fits(goodfit_ind,2,4);
    sigmay = lbub_fits(goodfit_ind,3,4);
    RFsize = sqrt(2*log(2))*geomean([sigmax sigmay],2);
    fprintf(', lbub_fits(ret)')
    
    % cellDists
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeTuneData.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_sizeTuneData.mat'] ' not found! Please remove from list\n'])
    end
    load(filename,'cellDists')
    fprintf(', cellDists')
    
    % sizeFitResults_SP
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_sizeFitResults_SP.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_sizeFitResults_SP.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'sizeFits')
    fprintf(', sizeFitResults_SP')
    
    % lbub_fits (size)
    filename = fullfile('\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind_size')
    fprintf(', lbub_fits(size) - done\n')
    
    switch expdata.area{iExp}
        case 'V1'
            RFcutoff = find(cellDists<10);
        case {'LM', 'AL'}
            RFcutoff = find(cellDists<15);
        case 'PM'
            RFcutoff = find(cellDists<20);
    end
    
    inds = intersect(goodfit_ind_size,RFcutoff);
    n = length(inds);
    
    %make maps
    retMap_El = NaN(size(mask_cell));
    retMap_Az = retMap_El;
    retMap_sigEl = retMap_El;
    retMap_sigAz = retMap_El;
    retMap_sigGeo = retMap_El;
    retMap_RFstim = retMap_El;
    retMap_prefSize = retMap_El;
    retMap_SI = retMap_El;
    % goodfit_ind = goodfit_ind(goodfit_ind_size);
    for i=1:length(goodfit_ind)
        ind = find(mask_cell == goodfit_ind(i));
        retMap_Az(ind) = lbub_fits_ret(goodfit_ind(i),4,4);
        retMap_El(ind) = lbub_fits_ret(goodfit_ind(i),5,4);
        retMap_sigGeo(ind) = RFsize(i);
        retMap_RFstim(ind) = cellDists(i);
        retMap_prefSize(ind) = sizeFits(i,end).prefSize;
        SI = sizeFits(i,end).suppInd;
        retMap_SI(ind) = SI;
        if SI>1; retMap_SI(ind)=1; elseif SI<0; retMap_SI(ind)=0; end
    end
    
    imAlpha=ones(size(retMap_El));
    imAlpha(isnan(retMap_El))=0; % set all unmasked pixels to alpha=0
    
    % plot maps
    figure(1);clf;
    colormap default
    %el
    subplot(2,3,1)
    imagesc(retMap_El,'AlphaData',imAlpha)
    title('goodfit cells: El')
    h = colorbar;
    ylabel(h,'El (deg)','Rotation',270.0,'VerticalAlignment','bottom')
    set(gca,'color',0*[1 1 1]);
    %az
    subplot(2,3,4)
    imagesc(retMap_Az,'AlphaData',imAlpha)
    title('goodfit cells: Az')
    h = colorbar;
    ylabel(h,'Az (deg)','Rotation',270.0,'VerticalAlignment','bottom')
    set(gca,'color',0*[1 1 1]);
    %sigGeo
    subplot(2,3,2)
    imagesc(retMap_sigGeo,'AlphaData',imAlpha)
    title('RF radius')
    h = colorbar;
    ylabel(h,'sigma_{Geo} (deg)','Rotation',270.0,'VerticalAlignment','bottom')
    set(gca,'color',0*[1 1 1]);
    %RFstim
    subplot(2,3,5)
    imagesc(retMap_RFstim,'AlphaData',imAlpha)
    title('RF-stim distance')
    h = colorbar;
    ylabel(h,'RF-SD (deg)','Rotation',270.0,'VerticalAlignment','bottom')
    set(gca,'color',0*[1 1 1]);
    %prefSize
    subplot(2,3,3)
    imagesc(retMap_prefSize,'AlphaData',imAlpha)
    title('prefSize')
    h = colorbar;
    ylabel(h,'pS (deg)','Rotation',270.0,'VerticalAlignment','bottom')
    set(gca,'color',0*[1 1 1]);
    %SI
    subplot(2,3,6)
    imagesc(retMap_SI,'AlphaData',imAlpha)
    title('SI')
    h = colorbar;
    ylabel(h,'SI','Rotation',270.0,'VerticalAlignment','bottom')
    set(gca,'color',0*[1 1 1]);
     
    suptitle(sprintf('Mouse: %s, date: %s, area: %s',mouse,date,expdata.area{iExp}))
    set(gcf, 'Position', [0,50,1600,800])
    close all
    %%
    figure(2)
    %imagesc(rot90(retMap_El),'AlphaData',rot90(imAlpha))
    imagesc(retMap_El,'AlphaData',(imAlpha))
    colormap jet
    colormap( flipud(oldmap) );
    h = colorbar;
    caxis([-15 15])
    %ylabel(h,'El (deg)','Rotation',270.0,'VerticalAlignment','bottom')
    h.Location = 'southoutside';
    ylabel(h,'El (deg)')
    set(gca,'Color',0*[1 1 1]);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    axis image
    set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off');
    %addpath('K:\Code\export_fig\altmany-export_fig-412662f')
    %export_fig('K:\ppts\_paper figs\fig1_cellEl2.pdf','-a2');
    print('N:\home\kevin\ppts\_paper figs\fig1 subfigs\fig1_cellEl.pdf','-dpdf')
    
    figure(3)
    %imagesc(rot90(retMap_Az),'AlphaData',rot90(imAlpha))
    imagesc((retMap_Az),'AlphaData',(imAlpha))
    colormap jet
    oldmap = colormap;
    colormap( flipud(oldmap) );
    h = colorbar;
    caxis([-10 20])
    %ylabel(h,'Az (deg)','Rotation',270.0,'VerticalAlignment','bottom')
    h.Location = 'southoutside';
    ylabel(h,'Az (deg)','HorizontalAlignment','center')
    set(gca,'color',0*[1 1 1]);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    axis image
    set(gcf,'Color',[1 1 1]); set(gca,'Color',[0 0 0]); set(gcf,'InvertHardCopy','off');
    print('K:\ppts\_paper figs\fig1_cellAz.pdf','-dpdf')
    
    if iExp == 39
        %pause
    end
end
%close all
