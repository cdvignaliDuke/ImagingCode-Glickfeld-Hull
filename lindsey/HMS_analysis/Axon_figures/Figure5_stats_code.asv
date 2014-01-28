fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 5;

P = 2;
matrix = 'SF5xTF5';
inj = 'LM';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

%average comparison
all_SFxTF_fit = zeros(5,5,3);
all_SFxTF_long = zeros(25,3);
for iArea = 1:3
    all_SFxTF_fit(:,:,iArea) = reshape(mean(Goodfits(iArea).plotfit,1),5,5)';  
    all_SFxTF_long(:,iArea) = reshape((all_SFxTF_fit(:,:,iArea)./max(max(all_SFxTF_fit(:,:,iArea),[],1),[],2))'*50,[25 1]);        
end

tuning_corr = triu2vec(corrcoef(all_SFxTF_long(:,1),all_SFxTF_long(:,3)));

%tuning curve comparisons
[h_PM_TF, p_PM_TF] = ttest(mean(Goodfits(1).plotfit(:,[1 2 6 7 11 12 16 17 21 22]),2), mean(Goodfits(1).plotfit(:,[4 5 9 10 14 15 19 20 24 25]),2));
[h_PM_SF, p_PM_SF] = ttest(mean(Goodfits(1).plotfit(:,1:10),2), mean(Goodfits(1).plotfit(:,16:25),2));
[h_PM_speed, p_PM_speed] = ttest(mean(Goodfits(1).plotfit(:,[1 2 3 6 7 11]),2), mean(Goodfits(1).plotfit(:,[15 19 20 23 24 25]),2));

[h_AL_TF, p_AL_TF] = ttest(mean(Goodfits(2).plotfit(:,[1 2 6 7 11 12 16 17 21 22]),2), mean(Goodfits(2).plotfit(:,[4 5 9 10 14 15 19 20 24 25]),2));
[h_AL_SF, p_AL_SF] = ttest(mean(Goodfits(2).plotfit(:,1:10),2), mean(Goodfits(2).plotfit(:,16:25),2));
[h_AL_speed, p_AL_speed] = ttest(mean(Goodfits(2).plotfit(:,[1 2 3 6 7 11]),2), mean(Goodfits(2).plotfit(:,[15 19 20 23 24 25]),2));

PM_SF_ratio = mean(mean(Goodfits(1).plotfit(:,1:10),2),1)./mean(mean(Goodfits(1).plotfit(:,16:25),2),1);
AL_SF_ratio = mean(mean(Goodfits(2).plotfit(:,16:25),2),1)./mean(mean(Goodfits(2).plotfit(:,1:10),2),1);

PM_TF_ratio = mean(mean(Goodfits(1).plotfit(:,[1 2 6 7 11 12 16 17 21 22]),2),1)./mean(mean(Goodfits(1).plotfit(:,[4 5 9 10 14 15 19 20 24 25]),2),1);
AL_TF_ratio = mean(mean(Goodfits(2).plotfit(:,[4 5 9 10 14 15 19 20 24 25]),2),1)./mean(mean(Goodfits(2).plotfit(:,[1 2 6 7 11 12 16 17 21 22]),2),1);

PM_speed_ratio = mean(mean(Goodfits(1).plotfit(:,[1 2 3 6 7 11]),2),1)./mean(mean(Goodfits(1).plotfit(:,[15 19 20 23 24 25]),2),1);
AL_speed_ratio = mean(mean(Goodfits(2).plotfit(:,[15 19 20 23 24 25]),2),1)./mean(mean(Goodfits(2).plotfit(:,[1 2 3 6 7 11]),2),1);

