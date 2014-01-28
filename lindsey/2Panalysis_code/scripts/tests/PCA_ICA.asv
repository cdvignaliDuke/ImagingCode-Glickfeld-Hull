%%
clear all
setpath;

global USER_PROFILE
USER_PROFILE = 'lindsey'; % 

protdir = dirs.images;
%%
base = 'G:\users\lindsey\regLG\';
date = '110509';
mouse = 'Y13';
userun = [1];
channel = '_green';

newdir = fullfile([date '_' mouse], ['run' num2str(userun)]);
expt = frGetExpt(newdir);

%% compute principal components           
newdir = '110902_Y18\run3';
frPrinComp(newdir);

%% load principal components
pcs = load(fullfile(expt.dirs.analrootpn,expt.filenames.pca_usv));
nt = size(pcs.v,1);

%% figure principal components
figure;
sm = stackFilter(pcs.U,1.5);

ax=[];
for pc = 1:16;
    ax(pc)=subplot(4,4,pc);
    imagesc(sm(:,:,pc));
    
    %imstretch(sm(:,:,pc),[.5 .99],1.5);
end;
colormap gray;

%% compute independent components
PCuse = [2:4 7:75];
mu = .5;
nIC = 32;
termtol = 1e-6;
maxrounds = 400;
mixedsig = pcs.v';
mixedfilters = pcs.U;
CovEvals = diag(pcs.s).^2;
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, ...
    mixedfilters, CovEvals, PCuse, mu, nIC,[],termtol,maxrounds);

% f0 = mean(reg,3);
dt = 1/frGetFrameRate;
tt = [0:nt-1]/frGetFrameRate;

%% plot independent filters
sel = 1:16;
cs = permute(ica_filters,[2,3,1]);
sm = stackFilter(cs,1.5);
figure;
ind = 1;
for ic = sel
    subplot(4,4,ind);
    imstretch(sm(:,:,ic),[.5 .99],1.5);
    ind = ind+1;
    text(.8,.1,num2str(ic),'fontsize',12,'color','w','fontweight','bold','unit','norm');
end;

%% plot independent signals
figure;
tcOffsetPlot(tt,ica_sig(1:4,:)',[], 'k-');


ica_sig_all = sum(ica_sig([1 2 3 4 5 6 7 8 9 12],:),1);
figure;plot(ica_sig_all)