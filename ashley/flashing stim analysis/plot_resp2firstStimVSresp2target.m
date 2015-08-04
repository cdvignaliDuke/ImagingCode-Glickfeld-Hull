edit FlashingStim_dataSortedByCycle_combineDatasetsSameType.m
% edit FlashingStim_dataSortedByCycle_combineDatasets.m
edit successTrials.m
edit resp2target.m

V_successAvg_target = V_successAvg;
AV_successAvg_target = AV_successAvg;

errbar_V_target = errbar_V;
errbar_AV_target = errbar_AV;

L = length(V_successAvg);

edit resp2firstStim.m

V_successAvg_base = V_successAvg;
AV_successAvg_base = AV_successAvg;
errbar_V_base = errbar_V;
errbar_AV_base = errbar_AV;

V_successAvg_targetNorm = V_successAvg_target - min(V_successAvg_target(30:34));
AV_successAvg_targetNorm = AV_successAvg_target - min(AV_successAvg_target(30:34));
V_successAvg_baseNorm = V_successAvg_base - min(V_successAvg_base(30:34));
AV_successAvg_baseNorm = AV_successAvg_base - min(AV_successAvg_base(30:34));

Lbase = length(V_successAvg_base);
%%

figure;
% shadedErrorBar([],V_successAvg_target,errbar_V_target,'g', 1);
% hold on
% shadedErrorBar([],AV_successAvg_target,errbar_AV_target,'m', 1);
% hold on
% shadedErrorBar([],V_successAvg_base,errbar_V_base,'g:', 1);
% hold on
% shadedErrorBar([],AV_successAvg_base,errbar_AV_base,'m:', 1);
plot(V_successAvg_target(30:end,:),'g');
hold on
plot(AV_successAvg_target(30:end,:),'m');
hold on
plot(V_successAvg_base(30:end,:),'g:');
hold on
plot(AV_successAvg_base(30:end,:),'m:');
hold on
hold on
vline(cycTime,'k:');

figure;
% shadedErrorBar([],V_successAvg_target,errbar_V_target,'g', 1);
% hold on
% shadedErrorBar([],AV_successAvg_target,errbar_AV_target,'m', 1);
% hold on
% shadedErrorBar([],V_successAvg_base,errbar_V_base,'g:', 1);
% hold on
% shadedErrorBar([],AV_successAvg_base,errbar_AV_base,'m:', 1);
plot(V_successAvg_targetNorm(30:end,:),'g');
hold on
plot(AV_successAvg_targetNorm(30:end,:),'m');
hold on
plot(V_successAvg_baseNorm(30:end,:),'g:');
hold on
plot(AV_successAvg_baseNorm(30:end,:),'m:');
hold on
hold on
vline(cycTime,'k:');