fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 3;

areas = strvcat('PM', 'LM', 'AL');
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

%distributions of TF, SF and speed
[h_PM_LM_TF, p_PM_LM_TF] = kstest2(log2(Goodfits(1).TF), log2(Goodfits(2).TF));
[h_PM_AL_TF, p_PM_AL_TF] = kstest2(log2(Goodfits(1).TF), log2(Goodfits(3).TF));
[h_AL_LM_TF, p_AL_LM_TF] = kstest2(log2(Goodfits(3).TF), log2(Goodfits(2).TF));

[h_PM_LM_SF, p_PM_LM_SF] = kstest2(log2(Goodfits(1).SF), log2(Goodfits(2).SF));
[h_PM_AL_SF, p_PM_AL_SF] = kstest2(log2(Goodfits(1).SF), log2(Goodfits(3).SF));
[h_AL_LM_SF, p_AL_LM_SF] = kstest2(log2(Goodfits(3).SF), log2(Goodfits(2).SF));

[h_PM_LM_speed, p_PM_LM_speed] = kstest2(log2(Goodfits(1).speed), log2(Goodfits(2).speed));
[h_PM_AL_speed, p_PM_AL_speed] = kstest2(log2(Goodfits(1).speed), log2(Goodfits(3).speed));
[h_AL_LM_speed, p_AL_LM_speed] = kstest2(log2(Goodfits(3).speed), log2(Goodfits(2).speed));

%difference in medians between mice
areas = strvcat('PM', 'LM', 'AL');
inj_list = {'V1' 'V1L5'};
matrix = 'SF5xTF5';
P = 2;
medians = [];
for iinj = 1:2
    inj = char(inj_list(iinj));
    if iinj == 1
        mouse_list = {'Y13' 'X32' 'DR7' 'DR9' 'AC39' 'AC42' 'AC44' 'AC45' 'Y18' 'Y26' 'M13' 'M14' 'M22' 'M31'};
    elseif iinj == 2
        mouse_list = {'VC1' 'CM94'};
    end
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
    load(fn_summary);
for iMouse = 1:length(mouse_list);
    mouse = mouse_list{iMouse};
    mouse_medians = [NaN NaN NaN];
    for iArea = 1:3
        image = areas(iArea,:);
        sum_base = 'G:\users\lindsey\analysisLG\experiments';
        list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
        load(list_fn);
        median_speed = [];
        nexp = all_fits(iArea).nexp;
        for iexp = 1:nexp
            exp_mouse = exp_list.mouse_mat{iexp};
            if length(exp_mouse)==length(mouse)
                if exp_mouse == mouse;
                    if all_fits(iArea).expt(iexp).n(2)>25
                        median_speed = [median_speed all_fits(iArea).expt(iexp).median_speed];
                    end
                end
            end
        end
        if length(median_speed)>0
            if length(median_speed)>1
                median_speed = mean(median_speed);
            end
            mouse_medians(:,iArea) = median_speed;
        end
    end
    medians = [medians; mouse_medians];
end
end

[h_PM_LM_medianspeed, p_PM_LM_mediansspeed] = ttest(medians(:,1), medians(:,2),[], 'left');
[h_PM_AL_mediansspeed, p_PM_AL_mediansspeed] = ttest(medians(:,1), medians(:,3),[], 'left');
[h_AL_LM_mediansspeed, p_AL_LM_mediansspeed] = ttest(medians(:,3), medians(:,2),[], 'right');

%comparison with cell bodies
[h_somas_boutons_PM, p_somas_boutons_PM] = kstest2(log2(Goodfits(1).speed), log2(Speed_PM_adj));
[h_somas_boutons_AL, p_somas_boutons_AL] = kstest2(log2(Goodfits(3).speed), log2(Speed_AL_adj));

x = floor(length(Goodfits(1).speed)./length(Speed_PM_adj));
y = rem(length(Goodfits(1).speed),length(Speed_PM_adj));
Speed_PM_mult = repmat(Speed_PM_adj, x,1);
Speed_PM_mult = [Speed_PM_mult; Speed_PM_adj(1:y,:)];
x = floor(length(Goodfits(3).speed)./length(Speed_AL_adj));
y = rem(length(Goodfits(3).speed),length(Speed_AL_adj));
Speed_AL_mult = repmat(Speed_AL_adj,x,1);
Speed_AL_mult = [Speed_AL_mult; Speed_AL_adj(1:y,:)];

p_kw = kruskalwallis([log2(Goodfits(1).speed) log2(Speed_PM_mult)]);
p_kw = kruskalwallis([log2(Goodfits(3).speed) log2(Speed_AL_mult)]);