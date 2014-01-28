%% reshuffle data
analdir = 'G:\users\lindsey\analysisLG\active mice\DR7\110505\';
seqfile = '110505_DR7_run6_max_Seqposition.mat';
load(fullfile(analdir,'analysis',seqfile));
[stims,blanks] = seq2epochs(Seqposition,18,18);

baseline = tcGetBaseline(ica_sig');
D = bsxfun(@minus,ica_sig',baseline);

[rr] = tcEpochAverage(D,stims);
[bb] =  tcEpochAverage(D,blanks);

sfs =[];

for iS = 1:length(Seqposition)
    sfs(iS)=Seqposition(iS).TFSFetc(2);
    cons(iS)=Seqposition(iS).TFSFetc(4);
end

%%
cs = permute(ica_filters,[2,3,1]);
sm = stackFilter(cs,1.5);

%% %%%%
figure; 
for ic = 1:6
   
r = nanmean(rr{ic});
b = nanmean(bb{ic});
re = nanstd(rr{ic});
be = nanstd(bb{ic});

usfs = unique(sfs);
nsfs = length(usfs);

for ind1 = 1:nsfs
    this = find(sfs==usfs(ind1) & cons>0);
       if ~isempty(this)
           arr1(ind1)=mean(r(this));
           arr2(ind1)=mean(b(this));
           arr3(ind1)=mean(re(this));
           arr4(ind1)=mean(be(this));
       end;
end;


axx=[];
axx(1)=subplot(3,3,ic);
errorbar(usfs, arr1, arr3,'ok-','markerfacecolor','k');
xlabel('Spatial Frequency (cpd)');
axis square
box off;
ylabel('Signal (a.u.)');
set(axx,'xlim',[0 .7])
end;


clf;
% ax=subplot(2,2,3);
% imagesc(arr1);colormap(flipud(gray));
% colorbar
% clim(ax,[0 .02]);
% axis(ax,'image');
% set(ax,'xtick',1:length(utfs),'xticklabel',utfs);
% set(ax,'ytick',1:length(usfs),'yticklabel',usfs);
% set(ax,'ydir','norm');
% xlabel('Temporal frequency (Hz)');
% ylabel('Spatial frequency (cpd)');



axx(2)=subplot(2,3,3);
errorbar(usfs,arr1(:,5),arr3(:,5),'sk-','markerfacecolor','k');
hold on;
errorbar(usfs,arr1(:,3),arr3(:,3),'sk-','markerfacecolor','w');
xlabel('Spatial frequency (cpd)');
axis square
legend('2 Hz','8 Hz');
box off;
set(axx,'ylim',[-.005 .03])
xlim([0.009,.64]);
set(gca,'xscale','log','xtick',lognums);
plot(xlim,r(end)*[1 1],'k:')

subplot(2,3,1);
imstretch(sm(:,:,ic),[.1 .99],1.5);
text(.8,.1,num2str(ic),'fontsize',12,'color','w','fontweight','bold','unit','norm');
title(sprintf('IC %i',ic))

stimulus = zeros(1,nt);
stimulus([stims{:}])=.95;
subplot(2,1,2);
hh(1)=area(tt2,stimulus*max(ica_sig(ic,:))*1.2)
hold on;
hh(2)=area(tt2,stimulus*min(ica_sig(ic,:))*1.2)
set(hh,'facecolor',.8*[1 1 1],'edgecolor','none');
plot(tt2,lowpass(ica_sig(ic,:)),'k');
hold on 
yl = ylim;
% plot(tt2,stimulus*yl(2),'r','linewidth',5);
xlim([0 400])
tt2 = [0:nt-1]/(frGetFrameRate/2);
xlabel('Time (s)');
box off;
% set(gcf,'paperposition',[0 0 5 5]);
% print(gcf,'-dpdf',fullfile(analdir,'figs',sprintf('axons_ic%i.pdf',ic)))

end;

