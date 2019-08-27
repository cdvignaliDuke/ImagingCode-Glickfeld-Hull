fn = 'img16_000_001';
nframes= 1000;
pn= 'Z:\Data\2P_imaging\150124img16';
pn_out = [pn '\analysis' fn '\'];

%load, baseline, register
cd = pn;
z= sbxread(fn,1,1000);        %Alterable line #frames
zmin = min(min(min(z,[],1),[],2),[],3);
zsub = z-zmin;
zavg = mean(zsub(:,:,1:50),3);    %rewrite this so that it selects a stable interval instead of just the first 50 frames
[zout, zreg] = stackRegister(zsub,zavg);
clear z
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
%save(fn_out,'-struct','pcs');
zregStruct = struct('zreg', zreg);     
fn_out2 = [pn_out [fn 'zreg.mat']]; 
%save(fn_out2, '-struct', 'zregStruct');   %zregStruct will save, it will just never leave the "busy" state. Check the folder to ensure it has saved then ctrl c


%visualize pca components
nt = size(pcs.v,1);
figure;
sm = stackFilter(pcs.U,1.5);
ax=[];
for pc = 1:25;                   % in order to visualize additional PCs simply alter the range (e.g. 26:50) Then subtract the appropriate amount from pc in the next line
    ax(pc)=subplot(5,5,pc);
    imagesc(sm(:,:,[pc]));
end;
colormap gray;

%compute independent components
PCuse = [1:50];
mu = 0;
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


%% TC amd ROI code
sel = [1:25];                           %changes numner of ICs graphed
cs = permute(ica_filters,[2,3,1]);
sm = stackFilter(cs,1.5);
figure;
ind = 1;
for ic = sel
    subplot(5,5,ind);                 %change here too
    imstretch(sm(:,:,ic),[.5 .99],1.5);
    ind = ind+1;
    text(.8,.1,num2str(ic),'fontsize',12,'color','w','fontweight','bold','unit','norm');
end;
clear mixedfilters pcs.u spatial_av spatial_norm mixedsig ica_filters ica_sig options defaultopts cs CovEvals ax ans zregStruct zout zsub zavg stack PCuse


%% threshold and extract TCs
sm_d = zeros(size(sm,1),size(sm,2),size(sel,2),3);      
sm_dsum = zeros(size(sm,1),size(sm,2),3);
%tc = zeros(size(zreg,3),size(sel,2));                    
b=1;
c= [.1, .5, .9];
c2 = NaN(size(sel));
for cc = 1:3:size(sel,2);
c2(1,[cc:(cc+2)]) = c;
end
tc_avg = zeros(size(zreg,3),length(sel));
b= 1;
for ic = sel
    [i,j] = find(sm(:,:,ic)>15e-13);   
    tc = zeros(size(i,1),size(zreg,3));
    for a = 1:size(i,1);
        sm_d(i(a),j(a),b,1) = b/size(sel,2); 
        sm_d(i(a),j(a),b,2) = c2(1,b);
        sm_d(i(a),j(a),b,3) = (size(sel,2)-b)/size(sel,2);
        tc(a,:) =(zreg(i(a),j(a),:));
   end
   tc_avg(:,b) = mean(tc,1)';
   for a =1:3;
   sm_dsum(:,:,a) = sm_dsum(:,:,a) + sm_d(:,:,b,a);    
   end
     b=b+1;
end

% for bb = 1:b-1;
%     figure;
%     imstretch(squeeze(sm_d(:,:,bb,1)));
% end

sm_dsum2 = sm_dsum;
for zz = 1:size(sm_dsum,3);
[xx,yy] = find(sm_dsum(:,:,zz)>1);
for aa = 1:size(xx,1);
    sm_dsum2(xx(aa),yy(aa),zz) = 1;
end 
end
figure; image(sm_dsum2); axis image; 
%figure; imstretch(sm_d(:,:,1,3),[.5 .99],1.5);
%figure; imstretch(sm_dsum(:,:),[.5 .99],1.5);
%figure; imagesc(sm_dsum); axis image; colormap(binary);
%figure; tcOffsetPlot2(tt,tc(:,[1:b-1]));

%% Plotting TCs

RGBmat = NaN(length(sel),3);
one = NaN(1,1);
two = NaN(1,1);
three = NaN(1,1);
for a = 1:size(RGBmat,1);
    one = a/size(sel,2);
    two = c2(1,a); 
    three = (size(sel,2)-a)/size(sel,2);
    RGBmat(a,:) = [one, two, three];
end
timeCourses = tc_avg(:,[1:b-1]);
delta = double(5*nanmean(nanstd(timeCourses)));
[nsamples,ncells]=size(timeCourses);
offsets = repmat(0:delta:(ncells-1)*delta,nsamples,1);
figure; hold; 
for ab = 1:size(RGBmat,1); 
    plot(tt,timeCourses(:,ab)+offsets(:,ab),'color', RGBmat(ab,:));
end
ylim([-delta (2.4*10^5)]);   %Hardcoded ylim xlim
xlim([0 100]);   

%% Saving
fn_out = [pn_out [fn '_TCs.mat']];
save(fn_out,'tc', 'sm_d', 'sm_dsum');


%%  Draw ROIs

% Alt method of selecting areas to draw
% I need to set all values above the threshold = 1
% set up a system that converts the edges of these zones to lines (1s)
% I have unthresholded ICs. Each on a diff figure
ROIdata = zeros(512,1250);   %NOTE should replace this with an image of avg fluorescence
sm_thresh = zeros(size(sm,1),size(sm,2),8);
a=1;
%sets everything above threshold to 1. creates a new matrix with all pixels
%either zero or one
for ic = sel;
    [i,j] = find(sm([2:511],[2:1249],ic)>7e-14);
    sm_thresh(i,j,a) = 1;
    a=a+1;
end

%makes the inner portion of the "ones" region = 0 so you end up with an
%outline
for w = 1:8;
    borders=[sm_thresh(i-1,j,w) + sm_thresh(i+1,j,w) + sm_thresh(i,j+1,w) + sm_thresh(i,j-1,w)];
    if borders == 4; 
        ROIdata(i,j) = 0;
    end
end
figure; imagesc(sm_thresh(:,:,3)); colormap(gray)
