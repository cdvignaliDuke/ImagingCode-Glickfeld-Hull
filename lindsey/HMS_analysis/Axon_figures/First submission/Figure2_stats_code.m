fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 2;
areas = strvcat('PM', 'LM', 'AL');
area_order = [2;3;1];
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
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

tuning_corr = zeros(1,3);
tuning_corr(:,1) = triu2vec(corrcoef(all_SFxTF_long(:,1),all_SFxTF_long(:,3)));
tuning_corr(:,2) = triu2vec(corrcoef(all_SFxTF_long(:,2),all_SFxTF_long(:,3)));
tuning_corr(:,3) = triu2vec(corrcoef(all_SFxTF_long(:,1),all_SFxTF_long(:,2)));

%mouse by mouse comparison
mouse_list = {'Y13' 'X32' 'DR9' 'AC39' 'AC44' 'AC45' 'Y26' 'M13' 'M14' 'M22' 'M31'};
mouse_tuning = zeros(25,11,3);
mouse_incl = zeros(length(mouse_list),3);
for iMouse = 1:length(mouse_list);
    mouse = mouse_list{iMouse};
    for iArea = 1:3
        image = areas(iArea,:);
        sum_base = 'G:\users\lindsey\analysisLG\experiments';
        list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
        load(list_fn);
        nexp = all_fits(iArea).nexp;
        plotfit = [];
        for iexp = 1:nexp
            exp_mouse = exp_list.mouse_mat{iexp};
            n = all_fits(iArea).expt(iexp).n(1);
            if length(exp_mouse)==length(mouse)
                if exp_mouse == mouse;
                    if all_fits(iArea).expt(iexp).n(2)>4
                        for iCell = 1:n;
                            if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1
                                plotfit = cat(3,plotfit, all_fits(iArea).expt(iexp).bouton(iCell).plotfit);
                            end
                        end
                    end
                end
            end
        end
        if size(plotfit,3)>4
        mouse_tuning(:,iMouse,iArea) = reshape(mean(plotfit,3)',25,1);
        mouse_incl(iMouse,iArea) = 1;
        end
    end
end

AL_PM_mouse_index = [1 4 5 6 7 8 9 10 11];
AL_LM_mouse_index = [5 6 7 8 9 10];
PM_LM_mouse_index = [3 5 6 7 8 9 10];

AL_PM_corr = zeros(1,length(AL_PM_mouse_index));
for iMouse = 1:length(AL_PM_mouse_index);
    AL_PM_corr(:,iMouse) = triu2vec(corrcoef(mouse_tuning(:,AL_PM_mouse_index(iMouse),1), mouse_tuning(:,AL_PM_mouse_index(iMouse),3)));
end
AL_LM_corr = zeros(1,length(AL_LM_mouse_index));
for iMouse = 1:length(AL_LM_mouse_index);
    AL_LM_corr(:,iMouse) = triu2vec(corrcoef(mouse_tuning(:,AL_LM_mouse_index(iMouse),2), mouse_tuning(:,AL_LM_mouse_index(iMouse),3)));
end
PM_LM_corr = zeros(1,length(PM_LM_mouse_index));
for iMouse = 1:length(PM_LM_mouse_index);
    PM_LM_corr(:,iMouse) = triu2vec(corrcoef(mouse_tuning(:,PM_LM_mouse_index(iMouse),2), mouse_tuning(:,PM_LM_mouse_index(iMouse),1)));
end

corr_avg = zeros(2,3);
corr_avg(1,1) = squeeze(mean(AL_PM_corr));
corr_avg(2,1) = squeeze(std(AL_PM_corr,[],2)./length(AL_PM_mouse_index));
corr_avg(1,2) = squeeze(mean(AL_LM_corr));
corr_avg(2,2) = squeeze(std(AL_LM_corr,[],2)./length(AL_LM_mouse_index));
corr_avg(1,3) = squeeze(mean(PM_LM_corr));
corr_avg(2,3) = squeeze(std(PM_LM_corr,[],2)./length(PM_LM_mouse_index));

%tuning curve comparisons
[h_PM_TF, p_PM_TF] = ttest(mean(Goodfits(1).plotfit(:,[1 2 6 7 11 12 16 17 21 22]),2), mean(Goodfits(1).plotfit(:,[4 5 9 10 14 15 19 20 24 25]),2));
[h_PM_SF, p_PM_SF] = ttest(mean(Goodfits(1).plotfit(:,1:10),2), mean(Goodfits(1).plotfit(:,16:25),2));
[h_PM_speed, p_PM_speed] = ttest(mean(Goodfits(1).plotfit(:,[1 2 3 6 7 11]),2), mean(Goodfits(1).plotfit(:,[15 19 20 23 24 25]),2));

[h_AL_TF, p_AL_TF] = ttest(mean(Goodfits(3).plotfit(:,[1 2 6 7 11 12 16 17 21 22]),2), mean(Goodfits(3).plotfit(:,[4 5 9 10 14 15 19 20 24 25]),2));
[h_AL_SF, p_AL_SF] = ttest(mean(Goodfits(3).plotfit(:,1:10),2), mean(Goodfits(3).plotfit(:,16:25),2));
[h_AL_speed, p_AL_speed] = ttest(mean(Goodfits(3).plotfit(:,[1 2 3 6 7 11]),2), mean(Goodfits(3).plotfit(:,[15 19 20 23 24 25]),2));

[h_LM_TF, p_LM_TF] = ttest(mean(Goodfits(2).plotfit(:,[1 2 6 7 11 12 16 17 21 22]),2), mean(Goodfits(2).plotfit(:,[4 5 9 10 14 15 19 20 24 25]),2));
[h_LM_SF, p_LM_SF] = ttest(mean(Goodfits(2).plotfit(:,1:10),2), mean(Goodfits(2).plotfit(:,16:25),2));
[h_LM_speed, p_LM_speed] = ttest(mean(Goodfits(2).plotfit(:,[1 2 3 6 7 11]),2), mean(Goodfits(2).plotfit(:,[15 19 20 23 24 25]),2));

PM_SF_ratio = mean(mean(Goodfits(1).plotfit(:,1:10),2),1)./mean(mean(Goodfits(1).plotfit(:,16:25),2),1);
AL_SF_ratio = mean(mean(Goodfits(3).plotfit(:,16:25),2),1)./mean(mean(Goodfits(3).plotfit(:,1:10),2),1);
LM_SF_ratio = mean(mean(Goodfits(2).plotfit(:,16:25),2),1)./mean(mean(Goodfits(2).plotfit(:,1:10),2),1);

PM_TF_ratio = mean(mean(Goodfits(1).plotfit(:,[1 2 6 7 11 12 16 17 21 22]),2),1)./mean(mean(Goodfits(1).plotfit(:,[4 5 9 10 14 15 19 20 24 25]),2),1);
AL_TF_ratio = mean(mean(Goodfits(3).plotfit(:,[4 5 9 10 14 15 19 20 24 25]),2),1)./mean(mean(Goodfits(3).plotfit(:,[1 2 6 7 11 12 16 17 21 22]),2),1);
LM_TF_ratio = mean(mean(Goodfits(2).plotfit(:,[4 5 9 10 14 15 19 20 24 25]),2),1)./mean(mean(Goodfits(2).plotfit(:,[1 2 6 7 11 12 16 17 21 22]),2),1);

PM_speed_ratio = mean(mean(Goodfits(1).plotfit(:,[1 2 3 6 7 11]),2),1)./mean(mean(Goodfits(1).plotfit(:,[15 19 20 23 24 25]),2),1);
AL_speed_ratio = mean(mean(Goodfits(3).plotfit(:,[15 19 20 23 24 25]),2),1)./mean(mean(Goodfits(3).plotfit(:,[1 2 3 6 7 11]),2),1);
LM_speed_ratio = mean(mean(Goodfits(2).plotfit(:,[15 19 20 23 24 25]),2),1)./mean(mean(Goodfits(2).plotfit(:,[1 2 3 6 7 11]),2),1);
