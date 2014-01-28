%%
clear all
setpath;

global USER_PROFILE
USER_PROFILE = 'lindsey'; % 

protdir = dirs.images;
%%
base = 'G:\users\lindsey\regLG\';
date = '110902';
mouse = 'Y18';
userun = [3:6];
channel = '_green';

%% register (single plane) and decimate
for iRun = 1:length(userun);
    newdir = fullfile([date '_' mouse], ['run' num2str(userun(iRun))]);
    expt = frGetExpt(newdir);
    frRegister(newdir,'Overwrite',true,'Oversampling',10,...
               'DoRecurse',false,'Engine','subpixel');
    outDir = fullfile('G:\users\lindsey\analysisLG\active mice\', mouse, date);
    if iRun == 1;            
        fn= fullfile(base, [date '_' mouse], ['run' num2str(userun(iRun))], ['run' num2str(userun(iRun)) channel]);
        stack = readtiff(fn);
        stack_down = stackGroupProject(stack, 12);
        av = mean(stack(:,:,1:5000),3);
        fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '_dec_reg.tif']);
        writetiff(stack_down, fn_out);
        else
        fn= fullfile(base, [date '_' mouse], ['run' num2str(userun(iRun))], ['run' num2str(userun(iRun)) channel]);
        stack = readtiff(fn);
        stack_down = stackGroupProject(stack, 12);
        [out, reg] = stackRegister(stack_down,av,10);
        fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '_dec_reg.tif']);
        writetiff(stack_down, fn_out);
    end
end
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
PCuse = [1 4:75];
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
tcOffsetPlot(tt,ica_sig(1:8,:)',[], 'k-');

ica_sig_down= zeros(nIC, 1320);
for time = 1:1320;
    ica_sig_down(:,time) = mean(ica_sig(:,1+(time-1)*10:10+(time-1)*10),2);
end
figure;
tt_dec = decimate(tt,10);
tcOffsetPlot(tt_dec, ica_sig_down(1:10,:)',[], 'k-');

figure;
ind =1;
for ic = [5 7 8 16 19]
    subplot(2,3,ind);
    imstretch(sm(:,:,ic),[.5 .99],1.5);
    ind = ind+1;
    text(.8,.1,num2str(ic),'fontsize',12,'color','w','fontweight','bold','unit','norm');
end;
