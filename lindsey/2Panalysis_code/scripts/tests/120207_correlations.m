P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120224';
mouse = 'X32';
date = '110512';
userun = [3:4];
nCond = 25;

base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse,date);

fn_resp = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp.mat']);
load(fn_resp);
fn_reps = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
load(fn_reps);
fn_lbub = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
load(fn_lbub);
fn_local = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
load(fn_local);

resp_dFoverF_shuf = zeros(size(resp_dFoverF));
start = 0;
for iCond = 1:nCond+1
    nRep = stim_reps(:,iCond);
    for iRep = 1:nRep-1
        resp_dFoverF_shuf(:,start+iRep+1)= resp_dFoverF(:,start+iRep);
    end
    resp_dFoverF_shuf(:,start+1) = resp_dFoverF(:,start+nRep);
    start = start+nRep;
end

r_real = zeros(size(resp_dFoverF,1));
r_shuf = zeros(size(resp_dFoverF,1));
for ipix1 = 1:size(resp_dFoverF,1)
    for ipix2 = 1:size(resp_dFoverF,1)
        r_real(ipix1,ipix2)= triu2vec(corrcoef(resp_dFoverF(ipix1,:),resp_dFoverF(ipix2,:)));
        r_shuf(ipix1,ipix2)= triu2vec(corrcoef(resp_dFoverF(ipix1,:),resp_dFoverF_shuf(ipix2,:)));
    end
end

r_norm = r_real-r_shuf;

r_real = zeros(length(goodfit_ind));
r_shuf = zeros(length(goodfit_ind));
for ipix1 = 1:length(goodfit_ind)
    for ipix2 = 1:length(goodfit_ind)
        r_real(ipix1,ipix2)= triu2vec(corrcoef(resp_dFoverF(goodfit_ind(ipix1),:),resp_dFoverF(goodfit_ind(ipix2),:)));
        r_shuf(ipix1,ipix2)= triu2vec(corrcoef(resp_dFoverF(goodfit_ind(ipix1),:),resp_dFoverF_shuf(goodfit_ind(ipix2),:)));
    end
end

r_norm = r_real-r_shuf;

figure; imagesq(r_real);
figure; imagesq(r_shuf);
figure; imagesq(r_norm);

[order, groups] = corrsort(r_real);
[order_shuf, groups_shuf] = corrsort(r_shuf);
[order_norm, groups_norm] = corrsort(r_norm);
figure;
subplot(2,3,1);
imagesq(r_real);
title('Real')
subplot(2,3,2);
imagesq(r_shuf);
title('Shuffled')
subplot(2,3,3);
imagesq(r_norm);
title('Real-Shuffled')
subplot(2,3,4);
imagesq(r_real(order,order));
subplot(2,3,5);
imagesq(r_shuf(order,order));
subplot(2,3,6);
imagesq(r_norm(order,order));

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_Y13_correlation_matrix.pdf']);
        print(gcf, '-dpdf', fn_out);

bouton_order = zeros(size(local_max));
bouton_order_shuf = zeros(size(local_max));
bouton_order_norm = zeros(size(local_max));
for ipix = 1:length(goodfit_ind);
    sub_y = [i(goodfit_ind(ipix))-1:i(goodfit_ind(ipix))+1]; 
    sub_x = [j(goodfit_ind(ipix))-1:j(goodfit_ind(ipix))+1];
    bouton_order(sub_y,sub_x) = find(order==ipix)+10;
    bouton_order_shuf(sub_y,sub_x) = find(order_shuf==ipix)+10;
    bouton_order_norm(sub_y,sub_x) = find(order_norm==ipix)+10;
end
figure; 
subplot(2,2,1)
imagesq(bouton_order);
title('Real')
subplot(2,2,2)
imagesq(bouton_order_shuf);
title('Shuf')
subplot(2,2,3)
imagesq(bouton_order_norm);
title('Norm')

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_X32_corrsort.pdf']);
        print(gcf, '-dpdf', fn_out);

%pca/kmeans        
SF_peak = [];
TF_peak = [];
SF_sigma = [];
TF_sigma = [];
xi = [];
n = all_fits(1).expt(2).n(1);        
for iCell = 1:n
    if all_fits(1).expt(2).bouton(iCell).goodfit == 1;
        SF_peak = [SF_peak all_fits(1).expt(2).bouton(iCell).SF_fit];
        TF_peak = [TF_peak all_fits(1).expt(2).bouton(iCell).TF_fit];
        SF_sigma = [SF_sigma all_fits(1).expt(2).bouton(iCell).sigma_SF];
        TF_sigma = [TF_sigma all_fits(1).expt(2).bouton(iCell).sigma_TF];
        xi = [xi all_fits(1).expt(2).bouton(iCell).xi_fit];
    end
end

params = [SF_peak' TF_peak' SF_sigma' TF_sigma' xi'];
idx = kmeans(params,40);
bouton_idx = zeros(size(local_max));
for ipix = 1:length(goodfit_ind);
    sub_y = [i(goodfit_ind(ipix))-1:i(goodfit_ind(ipix))+1]; 
    sub_x = [j(goodfit_ind(ipix))-1:j(goodfit_ind(ipix))+1];
    bouton_idx(sub_y,sub_x) = T(ipix)*2;
end
figure; imagesq(bouton_idx)


Y = pdist(params,'cityblock');
Z = linkage(Y,'average');
figure;
[H,T] = dendrogram(Z,'colorthreshold','default');
set(H,'LineWidth',2)

[coeff score] = princomp(params);
figure;
subplot(2,2,1)
scatter(score(:,2),score(:,1))
subplot(2,2,2)
scatter(score(:,2),score(:,3))
subplot(2,2,3)
scatter(score(:,2),score(:,4))
subplot(2,2,4)
scatter(score(:,2),score(:,5))

y = length(goodfit_ind);
x = jet(y);
cmap = x; (1:2:end,:);
figure;
for ipix = 1:length(goodfit_ind)
subplot(2,2,1)
scatter3(score(ipix,2),score(ipix,3), score(ipix,1), [], cmap(ipix,:));
title('X=2; Y=3; Z=1')
hold on
subplot(2,2,2)
scatter3(score(ipix,3),score(ipix,4), score(ipix,1), [], cmap(ipix,:))
title('X=3; Y=4; Z=1')
hold on
subplot(2,2,3)
scatter3(score(ipix,3),score(ipix,4), score(ipix,2), [], cmap(ipix,:))
title('X=3; Y=4; Z=2')
hold on
subplot(2,2,4)
scatter3(score(ipix,4),score(ipix,5), score(ipix,1), [], cmap(ipix,:))
title('X=4; Y=5; Z=1')
hold on
end
suptitle('Order color')
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_princomp_scatter_order.pdf']);
        print(gcf, '-dpdf', fn_out);
fn_out = fullfile(anal_base, [date '_' mouse '_run' num2str(userun) '_params.mat']);
save(fn_out, 'params', 'resp_dFoverF_fit', 'stim_reps');
      
resp_dFoverF_fit = resp_dFoverF(goodfit_ind,:);