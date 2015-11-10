ialign = 1;
awFSAVdatasets;
for iexp = 1:size(expt,2)
% iexp = 1;
divideupdatabyalignment
%%
%%

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);


%% find sets of cells
DirFolder = expt(iexp).dirtuning;
run('cellSets.m')
%%
frameratems = expt(iexp).frame_rate/1000;
cyctimems = cycTime/frameratems;
stimon_frames = input.nFramesOn;
stimoff_frames = input.nFramesOff;
stimontime = stimon_frames/frameratems;
stimofftime = stimoff_frames/frameratems;
trialresult = unique(trialOutcome);
responsedelay_frames = 3;
yminmax = [-0.1 1];

%%
cells = cellsSelect{1};
for i = 1:length(cells)
cell1 = cells(i);
CYC = 6;

cycInd = find(tCyclesOn == cycles(CYC));
successIx = find(ismember(cycInd,find(strcmp(trialOutcome,'success'))));
cycBlock2 = block2(cycInd);
tempdata = cycData{CYC};
tempdfoverf = cycDataDFoverF{CYC};

fSTD = std(tempdata(1:30,:,:),[],1);
f2xSTD = 2.*(std(tempdata(1:30,:,:),[],1));
fMean = mean(tempdata(1:30,:,:),1);
fMeanPlusF2xSTD = bsxfun(@plus,fMean,f2xSTD);
dataGT2xSTD = bsxfun(@gt,tempdata,fMeanPlusF2xSTD);
dfoverfGT2xSTD = zeros(size(tempdfoverf));
dfoverfGT2xSTD(dataGT2xSTD) = tempdfoverf(dataGT2xSTD);

%% how many stds (or 2x std) above the mean is the response at each frame for each cell
nStdGtBl = zeros(size(tempdata));
nStdGtBl = floor(bsxfun(@rdivide,bsxfun(@minus,tempdata,fMean),f2xSTD));
nStdGtBl(nStdGtBl<0) = 0;

eventsMeanAll_V = squeeze(mean(nanmean(nStdGtBl(:,:,intersect(find(cycBlock2 == 0),successIx)),3),2));
eventsMeanAll_AV = squeeze(mean(nanmean(nStdGtBl(:,:,intersect(find(cycBlock2 == 1),successIx)),3),2));
eventsSTDAll_V = squeeze(std(nanmean(nStdGtBl(:,:,intersect(find(cycBlock2 == 0),successIx)),3),[],2));
eventsSTDAll_AV = squeeze(std(nanmean(nStdGtBl(:,:,intersect(find(cycBlock2 == 1),successIx)),3),[],2));

%% find absolute number of events over time

nStdGtBl_diff = diff(nStdGtBl);
nStdGtBl_diff(nStdGtBl_diff<0) = 0;

eventsdiffMeanAll_V = [0; squeeze(mean(nanmean(nStdGtBl_diff(:,:,intersect(find(cycBlock2 == 0),successIx)),3),2))];
eventsdiffMeanAll_AV = [0; squeeze(mean(nanmean(nStdGtBl_diff(:,:,intersect(find(cycBlock2 == 1),successIx)),3),2))];
eventsdiffSTDAll_V = [0; squeeze(std(nanmean(nStdGtBl_diff(:,:,intersect(find(cycBlock2 == 0),successIx)),3),[],2))];
eventsdiffSTDAll_AV = [0; squeeze(std(nanmean(nStdGtBl_diff(:,:,intersect(find(cycBlock2 == 1),successIx)),3),[],2))];

%% bin events over trial time course
binsz = 3;
vlng = floor(length(eventsdiffMeanAll_V)/binsz);
eventsdiffmean_down_V = mean(reshape(eventsdiffMeanAll_V(1:vlng*binsz),binsz,vlng),1);
eventsdiffmean_down_AV = mean(reshape(eventsdiffMeanAll_AV(1:vlng*binsz),binsz,vlng),1);

%%
figure;
plot(eventsdiffmean_down_V,'g')
%%
figure;
xaxisMS = ([1:size(nStdGtBl,1)]-prepress_frames)/frameratems;

hold on
for i = 1:cycles(CYC)
    L = ((i-1)*cyctimems);
    cyclepatch = patch([L L+stimontime L+stimontime L],[yminmax(1) yminmax(1) yminmax(2) yminmax(2)],'k');
    set(cyclepatch,'FaceAlpha',0.15);
    set(cyclepatch,'EdgeColor','none');
    hold on
end
targetpatch = patch([(CYC*cyctimems) (CYC*cyctimems)+stimontime (CYC*cyctimems)+stimontime (CYC*cyctimems)],[yminmax(1) yminmax(1) yminmax(2) yminmax(2)],'k');
set(targetpatch,'FaceAlpha',0.45);
set(targetpatch,'EdgeColor','none');
hold on

xlim([xaxisMS(1) xaxisMS(end)]);
ylim(yminmax);


% errorbar(xaxisMS,eventsMeanAll_V,eventsSTDAll_V,'g')
plot(xaxisMS,eventsMeanAll_V,'g')
hold on
% errorbar(xaxisMS,eventsMeanAll_AV,eventsSTDAll_AV,'k')
plot(xaxisMS,eventsMeanAll_AV,'k')
hold on
plot(xaxisMS,eventsdiffMeanAll_V,'r.','MarkerSize',10)
hold on
plot(xaxisMS,eventsdiffMeanAll_AV,'b.','MarkerSize',10)

ylabel('events (2xSTD above BL)')
xlabel('frames')
title({[mouse '-' dateofexp '; avg events across cells and trials'];'gr=vis, blk=aud'})
try
    print([fnout '\eventdetectfigs\' mouse '-' dateofexp '_' num2str(CYC) 'cycles_nSTDs_avgallcells'], '-dpdf');
catch
    mkdir(fnout,'eventdetectfigs')
    print([fnout '\eventdetectfigs\' mouse '-' dateofexp '_' num2str(CYC) 'cycles_nSTDs_avgallcells'], '-dpdf');
end

%%
figure;
yminmax = [0 .5];
for i = 1:cycles(CYC)
    L = ((i-1)*cyctimems);
    cyclepatch = patch([L L+stimontime L+stimontime L],[yminmax(1) yminmax(1) yminmax(2) yminmax(2)],'k');
    set(cyclepatch,'FaceAlpha',0.15);
    set(cyclepatch,'EdgeColor','none');
    hold on
end
targetpatch = patch([(CYC*cyctimems) (CYC*cyctimems)+stimontime (CYC*cyctimems)+stimontime (CYC*cyctimems)],[yminmax(1) yminmax(1) yminmax(2) yminmax(2)],'k');
set(targetpatch,'FaceAlpha',0.45);
set(targetpatch,'EdgeColor','none');
hold on

xlim([xaxisMS(1) xaxisMS(end)]);
ylim(yminmax);

plot(xaxisMS(3:3:length(xaxisMS)),eventsdiffmean_down_V,'g.','MarkerSize',20)
hold on
plot(xaxisMS(3:3:length(xaxisMS)),eventsdiffmean_down_AV,'k.','MarkerSize',20)
ylabel('events (2xSTD above BL)')
xlabel('frames')
title({[mouse '-' dateofexp '; abs# events across cells and trials'];'gr=vis, blk=aud'})
try
    print([fnout '\eventdetectfigs\' mouse '-' dateofexp '_' num2str(CYC) 'cycles_absNevents_avgallcells'], '-dpdf');
catch
    mkdir(fnout,'eventdetectfigs')
    print([fnout '\eventdetectfigs\' mouse '-' dateofexp '_' num2str(CYC) 'cycles_absNevents_avgallcells'], '-dpdf');
end
end
%%
end
%%
% figure;
% subplot(1,2,1)
% tcOffsetPlot(squeeze(dataGT2xSTD(:,cell1,1:10)))
% % plot(squeeze(mean(dataGT2xSTD(:,cell1,:),3)))
% % ylim([-1 2])
% hold on
% vline(30,'k')
% for i = 1:CYC-1;
%     vline(30+(cycTime*i),'k:')
% end
% vline(30+(cycTime*CYC),'c')
% 
% 
% subplot(1,2,2)
% tcOffsetPlot(squeeze(tempdfoverf(:,cell1,1:10)))
% % plot(squeeze(mean(tempdfoverf(:,cell1,:),3)))
% % plot(squeeze(tempdfoverf(:,cell1,:)))
% hold on
% vline(30,'k')
% for i = 1:CYC-1;
%     vline(30+(cycTime*i),'k:')
% end
% vline(30+(cycTime*CYC),'c')
% 
% %%
% figure;
% % tcOffsetPlot(squeeze(dfoverfGT2xSTD(:,cell1,1:10)))
% % plot(squeeze(mean(dfoverfGT2xSTD(:,cell1,:),3)))
% tcOffsetPlot(mean(dfoverfGT2xSTD(:,drivencells,:),3))
% 
% hold on
% vline(30,'k')
% for i = 1:CYC-1;
%     vline(30+(cycTime*i),'k:')
% end
% vline(30+(cycTime*CYC),'c')