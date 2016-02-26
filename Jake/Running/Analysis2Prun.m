%{
Analysis2P for running movies. Adapted from Jake's Analysis2P.m.
Loads and performs PCA/ICA on 2P movie to create fluorescence time courses
for each detected component.
Minhal Ahmed 2/25/16
%}
clear;
pathname = ['S:\2PImaging\'];
moviename = ['160219img39\img39\img39_000_001'];
moviefile = [pathname moviename];
nframes = 6000;
%% Load movie, adjust baseline
movie2p = squeeze(sbxread(moviefile,1,nframes));
%Find minimum value in movie and subtract from original to set baseline to
%zero.
movie2pmin = min(min(min(movie2p,[],1),[],2),[],3);
movie2psub = movie2p - movie2pmin;
movie2pavg = mean(movie2psub(:,:,1:50),3); %Note: This arbitrarily selects the first 50 frames of the movie to serve as a "stable template." Consider improving.
[movie2pout, movie2preg] = stackRegister(movie2psub,movie2pavg);
writetiff(movie2p,'S:\2PImaging\160219img39\img39\img39_1')
%writetiff(movie2preg,'S:\2PImaging\160219img39\img39\img39_1reg')
clear movie2p
tiffmovie = ['S:\2PImaging\160219img39\img39\img39_1.tif'];

%% Prep for PCA

[mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(tiffmovie,[1 nframes], 300,1, [pathname '160219img39\img39'], []); 
[PCuse] = CellsortChoosePCs(tiffmovie, mixedfilters);
CellsortPlotPCspectrum(tiffmovie,CovEvals,PCuse);

%{
global stack
stack = single(movie2preg);
defaultopts = {'nComp', 300, 'BorderWidth',4};
options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);
[ny,nx,nt] = size(stack);
roi = imborder([ny,nx],options.BorderWidth,0); %creates a mask the size of the original movie w/ a border the size of borderWidth.
fprintf('Masking edges...');
stack = bsxfun(@times,stack,single(roi));
fprintf('Done\n');
%Compute thin svd using randomized algorithm
pcs = stackFastPCA(1,options.nComp);


%% Visualize PCA components

nt = size(pcs.v,1);
figure;
sm = stackFilter(pcs.U,1.5);
ax=[];
for pc = 1:25;
    ax(pc)=subplot(5,5,pc);
    imagesc(sm(:,:,[pc]));
end
colormap gray;
%}

%% Compute independent components
mu = 0.5;
nIC = 50;
termtol = 1e-6;
maxrounds = 400;
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, ...
    mixedfilters, CovEvals, PCuse, mu, nIC, [], termtol, maxrounds);
[ny,nx,nt]=size(movie2preg);
dt = 1/frGetFrameRate;tt = [0:nt-1]/frGetFrameRate;
CellsortICAplot('series',ica_filters,ica_sig,movie2pavg,[0 200],dt,tt,1,[1:50],[],[])


%% TC and ROI code
sel = [1:50];
cs = permute(ica_filters,[2,3,1]);
sm = stackFilter(cs,1.5);
figure;
ind = 1;
for ic = sel
    subplot(5,5,ind);
    imstretch(sm(:,:,ic),[.5 .99],1.5);
    ind = ind+1;
    text(.8,.1,num2str(ic),'fontsize',12,'color','w','fontweight','bold','unit','norm');
end
%clear mixedfilters pcs.u spatial_av spatial_norm mixedsig ica_filters ica_sig options defaultopts cs CovEvals ax ans zregStruct zout zsub zavg stack PCuse

%% Threshold and extract TCs
sm_d = zeros(size(sm,1),size(sm,2),size(sel,2),3);
sm_dsum = zeros(size(sm,1),size(sm,2),3);
b=1;
c = [.1, .5, .9]; c2 = NaN(size(sel));
for cc = 1:3:size(sel,2)
    c2(1,[cc:(cc+2)]) = c;
end
tc_avg = zeros(size(movie2preg,3),length(sel));
for ic=sel
    [i,j] = find(sm(:,:,ic)>0.1);
    tc = zeros(size(i,1),size(movie2preg,3));
    for a = 1:size(i,1);
        sm_d(i(a),j(a),b,1) = b/size(sel,2);
        sm_d(i(a),j(a),b,2) = c2(1,b);
        sm_d(i(a),j(a),b,3) = (size(sel,2)-b)/size(sel,2);
        tc(a,:) = (movie2preg(i(a),j(a),:));
    end
    tc_avg(:,b) = mean(tc,1)';
    for a = 1:3;
        sm_dsum(:,:,a) = sm_dsum(:,:,a) + sm_d(:,:,b,a);
    end
    b = b+1;
end

sm_dsum2 = sm_dsum;
for zz = 1:size(sm_dsum,3);
    [xx,yy] = find(sm_dsum(:,:,zz)>1);
    for aa = 1:size(xx,1);
        sm_dsum2(xx(aa),yy(aa),zz) = 1;
    end
end
figure;
image(sm_dsum2);
axis image;

%% Plotting TCs

RGBmat = NaN(length(sel),3);
one = NaN(1,1); two = NaN(1,1); three = NaN(1,1);
for a = 1:size(RGBmat,1);
    one = a/size(sel,2);
    two = c2(1,a);
    three = (size(sel,2)-a)/size(sel,2);
    RGBmat(a,:) = [one, two, three];
end
timeCourses = tc_avg(:,[1:b-1]);
delta = double(5*nanmean(nanstd(timeCourses)));
[nsamples, ncells] = size(timeCourses);
offsets = repmat(0:delta:(ncells-1)*delta,nsamples,1);
figure; hold;
for ab = 1:size(RGBmat,1);
    plot(tt,timeCourses(:,ab)+offsets(:,ab),'color',RGBmat(ab,:));
end
ylim([-delta (2.4*10^5)]);
xlim([0 100]);




