%% Concatenating training days into one matrix

elr = [];
elr = [earlyLeverReleaseTC_03 earlyLeverReleaseTC_04 earlyLeverReleaseTC_05 earlyLeverReleaseTC_06 earlyLeverReleaseTC_08 earlyLeverReleaseTC_09 earlyLeverReleaseTC_10 earlyLeverReleaseTC_11 earlyLeverReleaseTC_12 earlyLeverReleaseTC_13 earlyLeverReleaseTC_15 earlyLeverReleaseTC_17 earlyLeverReleaseTC_18 earlyLeverReleaseTC_19 earlyLeverReleaseTC_20];

stt_al = []; stt_cs = []; stt_cv = []; stt_lm = []; stt_rl = []; stt_v = [];
stt_al = [sucTrialTargetTC_AL_03 sucTrialTargetTC_AL_04 sucTrialTargetTC_AL_05 sucTrialTargetTC_AL_06 sucTrialTargetTC_AL_08 sucTrialTargetTC_AL_10 sucTrialTargetTC_AL_11 sucTrialTargetTC_AL_12 sucTrialTargetTC_AL_13 sucTrialTargetTC_AL_15 sucTrialTargetTC_AL_17 sucTrialTargetTC_AL_18 sucTrialTargetTC_AL_19 sucTrialTargetTC_AL_20];
stt_cs = [sucTrialTargetTC_CS_03 sucTrialTargetTC_CS_04 sucTrialTargetTC_CS_05 sucTrialTargetTC_CS_06 sucTrialTargetTC_CS_08 sucTrialTargetTC_CS_10 sucTrialTargetTC_CS_11 sucTrialTargetTC_CS_12 sucTrialTargetTC_CS_13 sucTrialTargetTC_CS_15 sucTrialTargetTC_CS_17 sucTrialTargetTC_CS_18 sucTrialTargetTC_CS_19 sucTrialTargetTC_CS_20];
stt_cv = [sucTrialTargetTC_CV_03 sucTrialTargetTC_CV_04 sucTrialTargetTC_CV_05 sucTrialTargetTC_CV_06 sucTrialTargetTC_CV_08 sucTrialTargetTC_CV_10 sucTrialTargetTC_CV_11 sucTrialTargetTC_CV_12 sucTrialTargetTC_CV_13 sucTrialTargetTC_CV_15 sucTrialTargetTC_CV_17 sucTrialTargetTC_CV_18 sucTrialTargetTC_CV_19 sucTrialTargetTC_CV_20];
stt_lm = [sucTrialTargetTC_LM_03 sucTrialTargetTC_LM_04 sucTrialTargetTC_LM_05 sucTrialTargetTC_LM_06 sucTrialTargetTC_LM_08 sucTrialTargetTC_LM_10 sucTrialTargetTC_LM_11 sucTrialTargetTC_LM_12 sucTrialTargetTC_LM_13 sucTrialTargetTC_LM_15 sucTrialTargetTC_LM_17 sucTrialTargetTC_LM_18 sucTrialTargetTC_LM_19 sucTrialTargetTC_LM_20];
stt_rl = [sucTrialTargetTC_RL_03 sucTrialTargetTC_RL_04 sucTrialTargetTC_RL_05 sucTrialTargetTC_RL_06 sucTrialTargetTC_RL_08 sucTrialTargetTC_RL_10 sucTrialTargetTC_RL_11 sucTrialTargetTC_RL_12 sucTrialTargetTC_RL_13 sucTrialTargetTC_RL_15 sucTrialTargetTC_RL_17 sucTrialTargetTC_RL_18 sucTrialTargetTC_RL_19 sucTrialTargetTC_RL_20];
stt_v = [sucTrialTargetTC_V_03 sucTrialTargetTC_V_04 sucTrialTargetTC_V_05 sucTrialTargetTC_V_06 sucTrialTargetTC_V_08 sucTrialTargetTC_V_10 sucTrialTargetTC_V_11 sucTrialTargetTC_V_12 sucTrialTargetTC_V_13 sucTrialTargetTC_V_15 sucTrialTargetTC_V_17 sucTrialTargetTC_V_18 sucTrialTargetTC_V_19 sucTrialTargetTC_V_20];

lp = [];
lp = [leverPressTC_V_03 leverPressTC_V_04 leverPressTC_V_05 leverPressTC_V_06 leverPressTC_V_08 leverPressTC_V_10 leverPressTC_V_11 leverPressTC_V_12 leverPressTC_V_13 leverPressTC_V_15 leverPressTC_V_17 leverPressTC_V_18 leverPressTC_V_19 leverPressTC_V_20];

mtt_al = []; mtt_cs = []; mtt_cv = []; mtt_lm = []; mtt_rl = []; mtt_v = [];
mtt_al = [missTrialTargetTC_AL_03 missTrialTargetTC_AL_04 missTrialTargetTC_AL_08 missTrialTargetTC_AL_11 missTrialTargetTC_AL_15 missTrialTargetTC_AL_17 missTrialTargetTC_AL_19];
mtt_cs = [missTrialTargetTC_CS_03 missTrialTargetTC_CS_04 missTrialTargetTC_CS_08 missTrialTargetTC_CS_11 missTrialTargetTC_CS_15 missTrialTargetTC_CS_17 missTrialTargetTC_CS_19];
mtt_cv = [missTrialTargetTC_CV_03 missTrialTargetTC_CV_04 missTrialTargetTC_CV_08 missTrialTargetTC_CV_11 missTrialTargetTC_CV_15 missTrialTargetTC_CV_17 missTrialTargetTC_CV_19];
mtt_lm = [missTrialTargetTC_LM_03 missTrialTargetTC_LM_04 missTrialTargetTC_LM_08 missTrialTargetTC_LM_11 missTrialTargetTC_LM_15 missTrialTargetTC_LM_17 missTrialTargetTC_LM_19];
mtt_rl = [missTrialTargetTC_RL_03 missTrialTargetTC_RL_04 missTrialTargetTC_RL_08 missTrialTargetTC_RL_11 missTrialTargetTC_RL_15 missTrialTargetTC_RL_17 missTrialTargetTC_RL_19];
mtt_v = [missTrialTargetTC_V_03 missTrialTargetTC_V_04 missTrialTargetTC_V_08 missTrialTargetTC_V_11 missTrialTargetTC_V_15 missTrialTargetTC_V_17 missTrialTargetTC_V_19];

slr_al = []; slr_cs = []; slr_cv = []; slr_lm = []; slr_rl = []; slr_v = [];
slr_al = [sucLeverReleaseTC_AL_03 sucLeverReleaseTC_AL_04 sucLeverReleaseTC_AL_05 sucLeverReleaseTC_AL_06 sucLeverReleaseTC_AL_08 sucLeverReleaseTC_AL_10 sucLeverReleaseTC_AL_11 sucLeverReleaseTC_AL_12 sucLeverReleaseTC_AL_13 sucLeverReleaseTC_AL_15 sucLeverReleaseTC_AL_17 sucLeverReleaseTC_AL_18 sucLeverReleaseTC_AL_19 sucLeverReleaseTC_AL_20];
slr_cs = [sucLeverReleaseTC_CS_03 sucLeverReleaseTC_CS_04 sucLeverReleaseTC_CS_05 sucLeverReleaseTC_CS_06 sucLeverReleaseTC_CS_08 sucLeverReleaseTC_CS_10 sucLeverReleaseTC_CS_11 sucLeverReleaseTC_CS_12 sucLeverReleaseTC_CS_13 sucLeverReleaseTC_CS_15 sucLeverReleaseTC_CS_17 sucLeverReleaseTC_CS_18 sucLeverReleaseTC_CS_19 sucLeverReleaseTC_CS_20];
slr_cv = [sucLeverReleaseTC_CV_03 sucLeverReleaseTC_CV_04 sucLeverReleaseTC_CV_05 sucLeverReleaseTC_CV_06 sucLeverReleaseTC_CV_08 sucLeverReleaseTC_CV_10 sucLeverReleaseTC_CV_11 sucLeverReleaseTC_CV_12 sucLeverReleaseTC_CV_13 sucLeverReleaseTC_CV_15 sucLeverReleaseTC_CV_17 sucLeverReleaseTC_CV_18 sucLeverReleaseTC_CV_19 sucLeverReleaseTC_CV_20];
slr_lm = [sucLeverReleaseTC_LM_03 sucLeverReleaseTC_LM_04 sucLeverReleaseTC_LM_05 sucLeverReleaseTC_LM_06 sucLeverReleaseTC_LM_08 sucLeverReleaseTC_LM_10 sucLeverReleaseTC_LM_11 sucLeverReleaseTC_LM_12 sucLeverReleaseTC_LM_13 sucLeverReleaseTC_LM_15 sucLeverReleaseTC_LM_17 sucLeverReleaseTC_LM_18 sucLeverReleaseTC_LM_19 sucLeverReleaseTC_LM_20];
slr_rl = [sucLeverReleaseTC_RL_03 sucLeverReleaseTC_RL_04 sucLeverReleaseTC_RL_05 sucLeverReleaseTC_RL_06 sucLeverReleaseTC_RL_08 sucLeverReleaseTC_RL_10 sucLeverReleaseTC_RL_11 sucLeverReleaseTC_RL_12 sucLeverReleaseTC_RL_13 sucLeverReleaseTC_RL_15 sucLeverReleaseTC_RL_17 sucLeverReleaseTC_RL_18 sucLeverReleaseTC_RL_19 sucLeverReleaseTC_RL_20];
slr_v = [sucLeverReleaseTC_V_03 sucLeverReleaseTC_V_04 sucLeverReleaseTC_V_05 sucLeverReleaseTC_V_06 sucLeverReleaseTC_V_08 sucLeverReleaseTC_V_10 sucLeverReleaseTC_V_11 sucLeverReleaseTC_V_12 sucLeverReleaseTC_V_13 sucLeverReleaseTC_V_15 sucLeverReleaseTC_V_17 sucLeverReleaseTC_V_18 sucLeverReleaseTC_V_19 sucLeverReleaseTC_V_20];

stt_al = []; stt_cs = []; stt_cv = []; stt_lm = []; stt_rl = []; stt_v = [];
stt_al = [sucTrialTargetTC_AL_03 sucTrialTargetTC_AL_04 sucTrialTargetTC_AL_05 sucTrialTargetTC_AL_20];
stt_cs = [sucTrialTargetTC_CS_03 sucTrialTargetTC_CS_04 sucTrialTargetTC_CS_05 sucTrialTargetTC_CS_20];
stt_cv = [sucTrialTargetTC_CV_03 sucTrialTargetTC_CV_04 sucTrialTargetTC_CV_05 sucTrialTargetTC_CV_20];
stt_lm = [sucTrialTargetTC_LM_03 sucTrialTargetTC_LM_04 sucTrialTargetTC_LM_05 sucTrialTargetTC_LM_20];
stt_rl = [sucTrialTargetTC_RL_03 sucTrialTargetTC_RL_04 sucTrialTargetTC_RL_05 sucTrialTargetTC_RL_20];
stt_v = [sucTrialTargetTC_V_03 sucTrialTargetTC_V_04 sucTrialTargetTC_V_05 sucTrialTargetTC_V_20];

    
%% Color gradation by training day

cmap = colormap('copper');
close
idx = round(linspace(1, size(cmap,1), 4));
cmap = cmap(idx,:);

figure; hold on;
for i_plt = 1:4
    plot(stt_cv(:,i_plt), '-', 'color', cmap(i_plt,:))
end

%% Subtracting minimum value to have similar baseline

for i = 1:14
    stt_new(:,i) = bsxfun(@minus, stt_al(:,i), min(stt_al(:,i)));
end
figure; plot(lp_new);
figure; plot(stt_V_633);

%% Plot mean F

for i  = 1:14
    F_TC_avg(i) = mean(F_TC{i});
end

%% Plot vertical line

y1 = get(gca, 'ylim');
x1 = 75;
plot([x1 x1],y1)