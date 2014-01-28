fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 6;

P = 2;
matrix = 'SF5xTF5';
inj = 'LM';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);

inj = 'LM';
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);
Goodfits_LM = Goodfits;

inj = 'V1';
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);
Goodfits_V1 = Goodfits;

%distributions of TF, SF and speed
[h_PM_AL_TF, p_PM_AL_TF] = kstest2(log2(Goodfits_LM(1).TF), log2(Goodfits_LM(2).TF));
[h_PM_AL_SF, p_PM_AL_SF] = kstest2(log2(Goodfits_LM(1).SF), log2(Goodfits_LM(2).SF));
[h_PM_AL_speed, p_PM_AL_speed] = kstest2(log2(Goodfits_LM(1).speed), log2(Goodfits_LM(2).speed));

%comparisons of V1 and LM
[h_LM_V1_PM_TF, p_LM_V1_PM_TF] = kstest2(log2(Goodfits_LM(1).TF), log2(Goodfits_V1(1).TF));
[h_LM_V1_PM_SF, p_LM_V1_PM_SF] = kstest2(log2(Goodfits_LM(1).SF), log2(Goodfits_V1(1).SF));
[h_LM_V1_PM_speed, p_LM_V1_PM_speed] = kstest2(log2(Goodfits_LM(1).speed), log2(Goodfits_V1(1).speed));

[h_LM_V1_AL_TF, p_LM_V1_AL_TF] = kstest2(log2(Goodfits_LM(2).TF), log2(Goodfits_V1(3).TF));
[h_LM_V1_AL_SF, p_LM_V1_AL_SF] = kstest2(log2(Goodfits_LM(2).SF), log2(Goodfits_V1(3).SF));
[h_LM_V1_AL_speed, p_LM_V1_AL_speed] = kstest2(log2(Goodfits_LM(2).speed), log2(Goodfits_V1(3).speed));
