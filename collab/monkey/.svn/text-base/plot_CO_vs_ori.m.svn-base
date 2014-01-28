function plot_CO_vs_ori (Oristats, COvalue)

alpha = 0.01;

% statistical significance
p_resp=[Oristats.p_value_resp];
p_sel=[Oristats.p_value_sel];

filter=find(p_resp<alpha);
filter2=find(p_sel<alpha);

% orientation selectivity
OSI=[Oristats.ori_vector_tune];
ori_width=[Oristats.ori_tuning_width];
ori_width(find(isnan(ori_width)))=90;
OI=[Oristats.OI];

% direction selectivity
DSI=[Oristats.dir_vector_tune];
Rmax_ori=[Oristats.R_best_ori];
dir_width=[Oristats.dir_tuning_width];
dir_width(find(isnan(dir_width)))=180;
DI=[Oristats.DI];
Rmax_dir=[Oristats.R_best_dir];

% preferred
best_ori=[Oristats.ori_vector_angle];
best_dir=[Oristats.best_dir_interp];


subplot(4,3,1);
plot(COvalue(filter),OSI(filter),'.')
title('CO vs OSI');

subplot(4,3,2);
plot(COvalue(filter),ori_width(filter),'.')
title('CO vs ori-width');

subplot(4,3,3);
plot(COvalue(filter),OI(filter),'.')
title('CO vs OI');

subplot(4,3,4);
plot(COvalue(filter),Rmax_ori(filter),'.')
title('CO vs Rmax-ori');

subplot(4,3,5);
plot(COvalue(filter),DSI(filter),'.')
title('CO vs DSI');

subplot(4,3,6);
plot(COvalue(filter),dir_width(filter),'.')
title('CO vs dir-width');

subplot(4,3,7);
plot(COvalue(filter),DI(filter),'.')
title('CO vs DI');

subplot(4,3,8);
plot(COvalue(filter),Rmax_dir(filter),'.')
title('CO vs Rmax-dir');

subplot(4,3,9);
plot(COvalue,-log10(p_resp),'.')
title('CO vs -log10(p-resp)');

subplot(4,3,10);
plot(COvalue,-log10(p_sel),'.')
title('CO vs -log10(p-sel)');

subplot(4,3,11);
plot(COvalue(filter2),best_ori(filter2),'.')
title('CO vs pref ori');

subplot(4,3,12);
plot(COvalue(filter2),best_dir(filter2),'.')
title('CO vs pref dir');
