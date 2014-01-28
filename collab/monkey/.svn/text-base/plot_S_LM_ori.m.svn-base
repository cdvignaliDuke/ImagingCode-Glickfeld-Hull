function plot_S_LM_ori (stats, Oristats)

contrast = [20.6, 8.4, 8.5, 8.8, 35.7, 9.7, 8.3, 7.3, 40, 80]/100;

Ncells=length(stats);
Ncolor=length(stats(1).dir_ratio_change);
ratio=reshape([stats.dir_ratio_change],Ncolor,Ncells);
nratio=ratio./repmat(contrast',1,Ncells);

% raw response (dF/F)
Sp=ratio(9,:);
Sm=ratio(10,:);
Mp=ratio(2,:);
Mm=ratio(6,:);
Lp=ratio(4,:);
Lm=ratio(8,:);
LpMp=ratio(1,:);
LmMm=ratio(5,:);
LpMm=ratio(7,:);
LmMp=ratio(3,:);

% normalized by contrast
nSp=nratio(9,:);
nSm=nratio(10,:);
nMp=nratio(2,:);
nMm=nratio(6,:);
nLp=nratio(4,:);
nLm=nratio(8,:);
nLpMp=nratio(1,:);
nLmMm=nratio(5,:);
nLpMm=nratio(7,:);
nLmMp=nratio(3,:);

% max response to cone stims
MaxS=max([Sp;Sm]);
nMaxS=max([nSp;nSm]);
MaxLM=max([Lp;Lm;Mp;Mm]);
nMaxLM=max([nLp;nLm;nMp;nMm]);
MaxLMS=max([Lp;Lm;Mp;Mm;Sp;Sm]);
nMaxLMS=max([nLp;nLm;nMp;nMm;nSp;nSm]);

% max respones to (L+M+,L-M-), or (L-M,M-L)
MaxsumLM=max([LpMp;LmMm]);
nMaxsumLM=max([nLpMp;nLmMm]);
MaxdifLM=max([LpMm;LmMp]);
nMaxdifLM=max([nLpMm;nLmMp]);

% orientation selectivity index
OSI=[Oristats.ori_vector_tune];

% show only responsive cells either color or ori
posi_color=statsFindPositiveResp(stats,0.01);
posi_ori=statsFindPositiveResp(Oristats,0.01);
filter=union(posi_color,posi_ori)

figure

subplot(5,4,1);
plot(MaxLMS(filter),OSI(filter),'.');
title('MaxLMS vs OSI');

subplot(5,4,2);
plot(nMaxLMS(filter),OSI(filter),'.');
title('MaxLMS(norm) vs OSI');

subplot(5,4,3);
plot(MaxS(filter),OSI(filter),'.');
title('MaxS vs OSI');

subplot(5,4,4);
plot(nMaxS(filter),OSI(filter),'.');
title('MaxS(norm) vs OSI');

subplot(5,4,5);
plot(MaxLM(filter),OSI(filter),'.');
title('MaxLM vs OSI');

subplot(5,4,6);
plot(nMaxLM(filter),OSI(filter),'.');
title('MaxLM(norm) vs OSI');

subplot(5,4,7);
plot(MaxdifLM(filter),OSI(filter),'.');
title('Max difLM vs OSI');

subplot(5,4,8);
plot(nMaxdifLM(filter),OSI(filter),'.');
title('Max difLM(norm) vs OSI');

subplot(5,4,9);
plot(MaxsumLM(filter),OSI(filter),'.');
title('Max sumLM vs OSI');

subplot(5,4,10);
plot(nMaxsumLM(filter),OSI(filter),'.');
title('Max sumLM(norm) vs OSI');

subplot(5,4,11);
plot(MaxS(filter),MaxLM(filter),'.');
title('MaxS vs MaxLM');

subplot(5,4,12);
plot(nMaxS(filter),nMaxLM(filter),'.');
title('MaxS(norm) vs MaxLM(norm)');

subplot(5,4,13);
plot(MaxLM(filter),MaxdifLM(filter),'.');
title('MaxLM vs Max difLM');

subplot(5,4,14);
plot(nMaxLM(filter),nMaxdifLM(filter),'.');
title('MaxLM(norm) vs Max difLM(norm)');

subplot(5,4,15);
plot(MaxLM(filter),MaxsumLM(filter),'.');
title('MaxLM vs Max sumLM');

subplot(5,4,16);
plot(nMaxLM(filter),nMaxsumLM(filter),'.');
title('MaxLM(norm) vs Max sumLM(norm)');

subplot(5,4,17);
plot(MaxS(filter),MaxdifLM(filter),'.');
title('MaxS vs Max difLM');

subplot(5,4,18);
plot(nMaxS(filter),nMaxdifLM(filter),'.');
title('MaxS(norm) vs Max difLM(norm)');

subplot(5,4,19);
plot(MaxS(filter),MaxsumLM(filter),'.');
title('MaxS vs Max sumLM');

subplot(5,4,20);
plot(nMaxS(filter),nMaxsumLM(filter),'.');
title('MaxS(norm) vs Max sumLM(norm)');





