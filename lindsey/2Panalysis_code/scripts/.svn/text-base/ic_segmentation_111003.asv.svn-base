n_ic = 6;
figure
for ic = 1:n_ic
    subplot(3,3,ic);
    imstretch(sm(:,:,ic),[.5 .99],1.5);
    text(.8,.1,num2str(ic),'fontsize',12,'color','w','fontweight','bold','unit','norm');
end;


[x y n] = size(sm);
thresh = zeros(x, y, n_ic);


figure;
for ic = 1:n_ic
m = max(max(sm(:,:,ic),[],2),[],1);
thresh(:,:,ic) = sm(:,:,ic)>0.1*m;
subplot(3,3, ic); imstretch(thresh(:,:,ic));
end



se1 = strel('disk',1);
ic_close = imclose(thresh,se1);
figure; 
for ic = 1:n_ic
subplot(3,3, ic); imstretch(ic_close(:,:,ic));
end
ic_open = imopen(ic_close,se1);
figure; 
for ic = 1:n_ic
subplot(3,3, ic); imstretch(ic_open(:,:,ic));
end

axon_mask = zeros(size(ic_open));
for ic = 1:n_ic
    axon_mask(:,:,ic) = bwlabel(ic_open(:,:,ic));
end
figure
for ic = 1:n_ic
    subplot(3,3, ic); imstretch(axon_mask(:,:,ic));
end

fn_stack = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_stack_uncorr.tif']);
stack_uncorr = readtiff(fn_stack);

avg_tc = mean(mean(stack_uncorr,2),1);
avg_f = mean(avg_tc,3);
avg_df = (avg_tc-avg_f)./avg_f;
figure;plot(squeeze(avg_df))

ind = find(avg_df>.05 | avg_df<-.05);
stack_uncorr_flash = stack_uncorr;
stack_uncorr_flash(:,:,ind) = [];
size(stack_uncorr_flash)
%clear stack_uncorr

max_segs = max(max(max(axon_mask,[],3),[],2),[],1);
timeCourses= zeros(size(stack_uncorr_flash,3), max_segs,n_ic);
for ic = 1:n_ic
   timeCourses(:,1:max(max(axon_mask(:,:,ic),[],2),[],1),ic) = stackGetTimeCourses(stack_uncorr_flash,axon_mask(:,:,ic));
end


figure;
r = zeros(max_segs,max_segs, ic);
for ic = 1:n_ic
    nRoi = max(max(axon_mask(:,:,ic),[],2),[],1);
    for iroi1 = 1:nRoi;
        for iroi2 = 1:nRoi;
            r(iroi1, iroi2,ic) = triu2vec(corrcoef(timeCourses(:,iroi1,ic),timeCourses(:,iroi2,ic)));
        end
    end
    subplot(3,3,ic);
    imstretch(r(:,:,ic));
    colormap hot;
    colorbar;
end

order = zeros(1,max_segs, n_ic);
groups = zeros(8,max_segs, n_ic);
figure;
for ic  = 1:n_ic
    [order(:,:,ic),groups(:,:,ic)] = corrsort(r(:,:,ic));
    subplot(3,3,ic)
    imagesc(r(order(:,:,ic), order(:,:,ic),ic));
    colormap hot;
end

%local correlations
b = 5;
r = zeros(x,y,max_segs,n_ic);
for ic = 1:n_ic
    nRoi = max(max(axon_mask(:,:,ic),[],2),[],1);
    for iroi = 1:nRoi;
        for iy = b+1:240-b
            fprintf('.');
            for ix = b+1:256-b
                sub = stack_uncorr_flash(iy-3:iy+3,ix-3:ix+3,:);
                sub_avg = mean(mean(sub,2),1);
                r(iy,ix,iroi,ic)= triu2vec(corrcoef(timeCourses(:,iroi,ic),sub_avg));
            end;
        end;
    end;
    figure;
    title(num2str(ic))
    ind = 1;
    nsubs = ceil(sqrt(nRoi));
    for iroi = 1:nRoi;
        subplot(nsubs,nsubs,ind); 
        imagesq(r(:,:,iroi,ic));
        colormap('hot');
        ind =ind+1;
    end
end



axon_reorder = zeros(x*y,n_ic);
for ic = 1:n_ic
    nRoi = max(max(axon_mask(:,:,ic),[],2),[],1);
    for imask = 1:nRoi;
        new_mask = find(axon_mask(:,:,ic) == imask);
        axon_reorder(new_mask,ic) = order(imask);
    end
end
axon_fig = reshape(axon_reorder,[x y n_ic]);
figure;
for ic = 1:n_ic
    subplot(3,3,ic)
    imstretch(axon_fig(:,:,ic));
    colormap hot
end

subplot(2,1,1);imagesc(lowpass(timeCourses(:,order)',[0 10]));colorbar
    