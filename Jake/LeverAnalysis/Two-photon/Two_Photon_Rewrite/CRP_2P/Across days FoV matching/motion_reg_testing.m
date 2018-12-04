%manually select two FoVs from each day. Try to grab the same neurons from
%an area that looks like a good match. Then register one to the other. 

load(['\\crash.dhe.duke.edu\data\home\jake\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\170524_000_img94\img_reg.mat']);
img_reg_d1 = img_reg;
load(['\\crash.dhe.duke.edu\data\home\jake\Analysis\Cue_reward_pairing_analysis\2P\PCA_ICA_recalibration\170529_000_img94\img_reg.mat']);
img_reg_PL = img_reg;

d1_dim1 = [5:174];
d1_dim2 = [300:572];%[150:647];
pl_dim1 = [1:170];
pl_dim2 = [306:578];%[156:653];

%img_ref_d1 = max(img_reg_d1(d1_dim1,d1_dim2,[1:2000]),[],3);
img_ref_d1 = mean(img_reg_d1(d1_dim1,d1_dim2,[1:2000]),3);
PL_frame_num = size(img_reg_d1,3);
[~, reg_match_PL] = stackRegister(img_reg_PL(pl_dim1,pl_dim2, [1:round(PL_frame_num/8)]), img_ref_d1);
%[~, reg_match_PL] = stackRegister(   mean(img_reg_PL(pl_dim1,pl_dim2, [1:2000]),3)  , img_ref_d1);

figure; 
subplot(2,1,1);
imagesc(img_ref_d1);
title('day1 max');

subplot(2,1,2);
%imagesc(max(reg_match_PL(:,:,[1:2000]),[],3));
imagesc(mean(reg_match_PL(:,:,[1:2000]),3));
%imagesc(mean(reg_match_PL,3));
title('PL max');
suptitle('img94 stack register');
%aa = corr2(img_ref_d1, max(reg_match_PL(:,:,[401:800]),[],3))
aa = corr2(img_ref_d1, mean(reg_match_PL(:,:,[401:800]),3))
%aa = corr2(img_ref_d1, reg_match_PL)
 
% PL_max_proj = max(img_reg_PL(  pl_dim1, pl_dim2,[1:2000]),[],3); 
% PL_max_proj =  imhistmatch(    PL_max_proj,   img_ref_d1        );
% [~,PL_reg_dem] = imregdemons(PL_max_proj, img_ref_d1);
%PL_mean_proj = mean(img_reg_PL(  pl_dim1, pl_dim2,[1:2000]),3); 
PL_mean_proj = mean(reg_match_PL(  :, :,[1:2000]),3); 
%PL_mean_proj = reg_match_PL; 
%PL_mean_proj =  imhistmatch(    PL_mean_proj,   img_ref_d1        );
[~,PL_reg_dem] = imregdemons(PL_mean_proj, img_ref_d1);

figure; 
subplot(2,1,1);
imagesc(img_ref_d1);
title('day1 max');

subplot(2,1,2);
imagesc(PL_reg_dem);
title('PL max');
 aa = corr2(img_ref_d1, PL_reg_dem)
 
 figure; imshowpair(img_ref_d1, PL_reg_dem);
% 
% subplot(3,1,3);
% imagesc(PL_reg_dem(:,:,2));
% title('PL max');
% aa = corr2(img_ref_d1, PL_reg_dem(:,:,2))
suptitle('img94 imregdemons')
%% Candidates for matching neurons across days and the coordinates of the FoV which is a best fit 


%img90
%D1 X[131 573]  Y[20:202]
%PL X[135:577]  Y[20:202]

%img91
%D1  x[227:517] Y[68 225]
%PL  X[257:547] Y[84:241]

%img92
%D1 X [83:332]  Y[25:252]
%PL X [76:325]  Y[35:262]

%img93
%D1 x[238 512] y[1 135]
%pl x[252 526] y[1 135]

%img89
%D1 x[183 526] y[0 132]
%PL x[220 563] y[0 132]

%img94
%D1 x[150:647] y[5:174]
%PL x[156:653] y[1:170]

%img036
%D1 x[208 494] y[155 264]
%PL x[215 501] y[155 264]

%img044
%D1 x[53 393] y[9 158]
%PL x[63 403] y[9 158]

%img050
%D1 x[1 386] y[1 264]
%PL x[1 386] y[1 264]

%img067
%D1 x[332 565] y[12 167]
%PL x[353 586] y[9 164]

%img070  all pixels
%D1 x[] y[]
%PL x[] y[]

%img077 None

%img081
%D1 X [100:381] Y [50:250]
%PL X [35:316]  Y [50:250]

%img084
%D1 x[209 506] y[48 197]
%PL x[204 511] y[34 183]

%img085  None


