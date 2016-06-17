%{
Analysis2P for running movies. Adapted from Jake's Analysis2P.m.
Loads and performs PCA/ICA on 2P movie to create fluorescence time courses
for each detected component.
Minhal Ahmed 2/25/16
%}
clear;
pathname = ['S:\2PImaging\160228_img40\img40\'];
moviename = ['img40_000_000'];
moviefile = [pathname moviename '.sbx'];
nframes = 30849;
%% Load movie, adjust baseline
%tiffmovie = [pathname 'img49.tif'];
load([pathname 'movies.mat'],'movie2preg','movie2pout');
load([pathname 'rundata.mat']);
frameIntMs = 33;
movie2pavg = mean(movie2preg(:,:,1:50),3);
global rt
rt = movie2preg;
clear rttiff;
%% Prep for PCA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prep for pca
global stack
stack = single(movie2preg);
[ny,nx,nt]=size(stack);
% compute thin svd using randomized algorithm
pcs = stackFastPCA(1,300);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculate principle components using CellsortPCA2 and choose visualized
%components.
[mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA2(moviefile,[1 nframes], 300,1, [pathname], []); 
[PCuse] = CellsortChoosePCs(movie2pavg, mixedfilters);

%Allow user input to remove bad PCs.
badpcs = input('Enter the values of unwanted PCs (ex. [7 8 14]): ');
for ii = 1:length(badpcs)
    PCuse = PCuse(PCuse~=badpcs(ii));
end

%figure;
%CellsortPlotPCspectrum(tiffmovie,CovEvals,PCuse);


%% Compute independent components
mu = 0.25;
nIC = input('Enter the desired number of ICs (must be < PCs): ');
termtol = 1e-6;
maxrounds = 400;
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, ...
    mixedfilters, CovEvals, PCuse, mu, nIC, [], termtol, maxrounds);

save([pathname 'newPC_1'], 'PCuse' ,'mixedsig' ,'mixedfilters', 'CovEvals', 'covtrace' ,'movm', 'movtm','ica_sig', 'ica_filters', 'ica_A', 'numiter')


[ny,nx,nt]=size(movie2preg);
dt = 1/frGetFrameRate;tt = [0:nt-1]/frGetFrameRate;

%{
thresh = 2.5;
[spmat, spt, spc] = CellsortFindspikes(ica_sig,thresh,dt,0,1);
figure;
CellsortICAplot('series',ica_filters,ica_sig,movie2pavg,[0 200],dt,tt,1,[1:20],spt, spc)
figure;
CellsortICAplot('series',ica_filters,ica_sig,movie2pavg,[0 200],dt,tt,3,[1:20],spt, spc)
figure;
%}

%Calculate and plot running data
B = [1:length(dataStruct.locomotionMatFrames(3:end))];
B = B*(frameIntMs/1000);
runvFfig = figure;
%Set/calculate variables with proper units
revoLength = 30;revoPulse = 128;    %revoLength = length of 1 revolution (cm). revoPulse = # of pulses/revolution.
runSpeed = (1000/frameIntMs)*(revoLength/revoPulse)*dataStruct.locomotionMatFrames(3:end); %runSpeed=(frames/sec)*(cm/pulse)*(pulses/frame)
runyes = runSpeed;
runyes(abs(runSpeed)>=4)=median(runSpeed(runSpeed>4));runyes(abs(runSpeed)<4)=0;
plot(B,runSpeed,B,runyes,'r')
title('Running Speed vs. Fluorescence')
ylabel('Running Speed (cm/s)')
xlabel('Time (s)')

sel = [1:nIC];
cs = permute(ica_filters,[2,3,1]);
sm = stackFilter(cs,1.5);

%% Select ROIs from ICs
mask_cell = zeros(size(sm));
sm_logical = zeros(264,796);

for ic = sel
    sm_logical((sm(:,:,ic)>mean([max(prctile(sm(:,:,ic),96,1)) max(prctile(sm(:,:,ic),96,2))])))=1;
    sm_logical((sm(:,:,ic)<=mean([max(prctile(sm(:,:,ic),96,1)) max(prctile(sm(:,:,ic),96,2))])))=0;
    sm_logical = logical(sm_logical);
    mask_cell(:,:,ic) = bwlabel(sm_logical);
end

 %consolidates all ROIs within IC into single ROI
    thresh = 0.94; %correlation threshold for calling two dendrites one thing
    sz = size(movie2preg);
    mask_cell_temp = zeros(sz(1)*sz(2), nIC);
    for ic = sel
        if length(unique(reshape(mask_cell(:,:,ic),[1 sz(1)*sz(2)])))>2
            data_tc_temp = stackGetTimeCourses(movie2preg(:,:,5000:8000),mask_cell(:,:,ic));
            data_corr_temp = corrcoef(data_tc_temp);
            ind_rem = 1:length(unique(reshape(mask_cell(:,:,ic),[1 sz(1)*sz(2)])))-1;
            for i = 1:length(unique(reshape(mask_cell(:,:,ic),[1 sz(1)*sz(2)])))-1
                ind = ind_rem(find(data_corr_temp(min(ind_rem,[],2),ind_rem)>thresh));
                if length(ind)>1
                    for ii = ind
                        if i == 1
                            mask_cell_temp(find(mask_cell(:,:,ic)== ii),ic) = 1;
                        else
                            mask_cell_temp(find(mask_cell(:,:,ic)== ii),nIC) = 1;
                        end
                    end
                else
                    if i == 1
                        mask_cell_temp(find(mask_cell(:,:,ic)== i),ic) = 1;
                    else
                        mask_cell_temp(find(mask_cell(:,:,ic)== ind),nIC) = 1;
                    end  
                end
                ind_rem = ind_rem(~ismember(ind_rem,ind));
                if ~isempty(ind_rem)
                    cat(3, mask_cell_temp, zeros(size(mask_cell_temp(:,:,1))));
                    nIC = nIC+1;
                else
                    break
                end
            end
        else
            mask_cell_temp(find(mask_cell(:,:,ic)),ic) = 1;
        end
    end
    
 %get preliminary timecourses for segregating and grouping ROIs
    data_tc = zeros(size(movie2preg(:,:,5000:8000),3), nIC);
    for ic = sel;
        if sum(mask_cell_temp(:,ic),1)>0
            data_tc(:,ic) = stackGetTimeCourses(movie2preg(:,:,5000:8000), reshape(mask_cell_temp(:,ic), [sz(1) sz(2)]));
        end
    end
    data_corr = corrcoef(data_tc);
    figure; imagesc(data_corr)
    
 %consolidate timecourses that are highly correlated
    [i j] = find(and(data_corr>thresh,data_corr<1));
    n = size(i,1);
    if n>1
        for ii = 1:(n/2)
            ind = find(mask_cell_temp(:,i(ii)));
            mask_cell_temp(ind,j(ii)) = 1;
            mask_cell_temp(ind,i(ii)) = 0;
        end
    end 
    
 %finds overlapping pixels of ROIs and based on correlations decides whether
    %to group them or to split them- if splitting, then overlapping pixels are
    %eliminated from both ROIs

    mask_overlap = zeros(1,sz(1)*sz(2));
    mask_all = zeros(1,sz(1)*sz(2));
    count = 0;
    for ic = 1:nIC
        ind_new = find(mask_cell_temp(:,ic))';
        if length(ind_new)>1
            ind_old = find(mask_all);
            overlap = ismember(ind_old,ind_new);
            ind_both = find(overlap);
            if length(ind_both)>1
                ic_match = unique(mask_all(ind_old(ind_both)));
                for im = 1:length(ic_match)
                    if data_corr(ic, ic_match(im))> 0.95
                        count = count+1;
                        mask_all(ind_new) = ic_match(im);
                    else
                        mask_all(ind_new) = ic;
                        mask_all(ind_old(ind_both)) = 0;
                        mask_overlap(ind_old(ind_both)) = 1;
                    end
                end
            else
                 mask_all(ind_new) = ic;
            end
        end
    end
    
 % removes ICs smaller than 200 pixels, renumbers ROIs so in continuous ascending order
    start = 1;
    mask_final = zeros(size(mask_all));
    for ic = 1:max(mask_all,[],2)
        ind = find(mask_all==ic);
        if length(ind)<300
             mask_overlap(find(mask_all==ic)) = 1;
             mask_all(ind) = 0;
        end
        ind = find(mask_all==ic);
        if length(ind)>0
            mask_final(ind)=start;
            start= start+1;
        end
    end

    data_tc = stackGetTimeCourses(movie2preg, reshape(mask_final,[sz(1) sz(2)]));
    
%% Threshold and extract TCs
nmask = max(max(mask_final));
%nmask = size(sel,2);
sm = zeros(size(sm,1),size(sm,2),nmask);
for ii=1:nmask
    sm_placeholder = zeros(size(sm,1),size(sm,2));
    sm_placeholder(mask_final==ii) = 1;
    sm(:,:,ii)= sm_placeholder;
end
        
sm_d = zeros(size(sm,1),size(sm,2),nmask,3);
sm_dsum = zeros(size(sm,1),size(sm,2),3);
b=1;
c = [.1, .5, .9]; c2 = NaN(1, nmask);
for cc = 1:3:nmask
    c2(1,[cc:(cc+2)]) = c;
end
tc_avg = zeros(size(movie2preg,3),nmask);
for ic=1:nmask
    [i,j] = find(sm(:,:,ic)==1);                             %>mean([max(prctile(sm(:,:,ic),96,1)) max(prctile(sm(:,:,ic),96,2))]));
    tc = zeros(size(i,1),size(movie2preg,3));
    for a = 1:size(i,1);
        sm_d(i(a),j(a),b,1) = b/nmask;
        sm_d(i(a),j(a),b,2) = c2(1,b);
        sm_d(i(a),j(a),b,3) = (nmask-b)/nmask;
        tc(a,:) = (movie2preg(i(a),j(a),:)); 
    end
    sorttc = sort(tc,2);                                    %Finds a baseline F for each mask by taking the average of the 200 lowest F values for each pixel in the mask.
    tcbaseF = mean(mean(tc(1:size(i,1),1:200)));                
    tc_avg(:,b) = mean(tc,1)';
    tc_avg(:,b) = (tc_avg(:,b)./tcbaseF)-1;
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

RGBmat = NaN(nmask,3);
one = NaN(1,1); two = NaN(1,1); three = NaN(1,1);
for a = 1:size(RGBmat,1);
    one = a/nmask;
    two = c2(1,a);
    three = (nmask-a)/nmask;
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
ylim([-delta 25]);
xlim([0 100]);

