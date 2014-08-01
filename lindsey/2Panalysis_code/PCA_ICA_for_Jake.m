fn = 'img11_000_000';
nframes= 1500;
pn= 'D:\Jake_temp\';
pn_out = 'D:\Jake_temp\analysis\';

%load, baseline, register
z= squeeze(readsbx([pn fn],0,1500));
zmin = min(min(min(z,[],1),[],2),[],3);
zsub = z-zmin;
zavg = mean(zsub(:,:,1:50),3);
[zout, zreg] = stackRegister(zsub,zavg);

%prep for pca
global stack
stack = single(zreg);
defaultopts = {'nComp',300,'BorderWidth',4};
options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);
[ny,nx,nt]=size(stack);
roi = imborder([ny,nx],options.BorderWidth,0); 
fprintf('Masking edges... ');
stack= bsxfun(@times,stack,single(roi));
fprintf('Done\n');
% compute thin svd using randomized algorithm
pcs = stackFastPCA(1,options.nComp);
% save extracted components 
fprintf('Saving principal components\n');
fn_out = [pn_out [fn '_pca_usv.mat']];
save(fn_out,'-struct','pcs');

%visualize pca components
nt = size(pcs.v,1);
figure;
sm = stackFilter(pcs.U,1.5);
ax=[];
for pc = 1:16;
    ax(pc)=subplot(4,4,pc);
    imagesc(sm(:,:,pc));
end;
colormap gray;

%compute independent components
PCuse = [1:50];
mu = .5;
nIC = 32;
termtol = 1e-6;
maxrounds = 400;
mixedsig = pcs.v';
mixedfilters = pcs.U;
CovEvals = diag(pcs.s).^2;
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, ...
    mixedfilters, CovEvals, PCuse, mu, nIC,[],termtol,maxrounds);

dt = 1/frGetFrameRate;
tt = [0:nt-1]/frGetFrameRate;

%plot independent filters
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

%threshold ICs and extract TCs
sm_d = zeros(size(sm));
sm_dsum = zeros(size(sm,1),size(sm,2));
tc = zeros(size(zreg,3),32);
for ic = 1:32
    [i,j] = find(sm(:,:,ic)>7e-13);
    sm_d(i,j,ic) = 1;
    sm_dsum = sm_dsum+ (sm_d(:,:,ic).*ic);
    tc(:,ic) = squeeze(mean(mean(zreg(i,j,:),1),2));
end
figure; imagesc(sm_dsum)
figure; tcOffsetPlot(tt,tc);

fn_out = [pn_out [fn '_TCs.mat']];
save(fn_out,'tc', 'sm_d', 'sm_dsum');
