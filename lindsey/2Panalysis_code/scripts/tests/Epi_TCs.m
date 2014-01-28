figure;
subplot(2,3,1)
errorbar(TC_avg(:,5,1,1),TC_sem(:,5,1,1),'c')
hold on
errorbar(TC_avg(:,5,13,1),TC_sem(:,5,13,1),'g')
hold on
errorbar(TC_avg(:,5,25,1),TC_sem(:,5,25,1),'r')
ylim([-0.25 1.5])
title ('PM')
subplot(2,3,2)
errorbar(TC_avg(:,1,1,1),TC_sem(:,1,1,1),'c')
hold on
errorbar(TC_avg(:,1,13,1),TC_sem(:,1,13,1),'g')
hold on
errorbar(TC_avg(:,1,25,1),TC_sem(:,1,25,1),'r')
ylim([-0.25 1.5])
title ('AL')
subplot(2,3,3)
errorbar(TC_avg(:,2,1,1),TC_sem(:,2,1,1),'c')
hold on
errorbar(TC_avg(:,2,13,1),TC_sem(:,2,13,1),'g')
hold on
errorbar(TC_avg(:,2,25,1),TC_sem(:,2,25,1),'r')
ylim([-0.25 1.5])
title ('RL')

subplot(2,3,4)
errorbar(TC_avg(:,3,1),TC_sem(:,3,1),'c')
hold on
errorbar(TC_avg(:,3,9),TC_sem(:,3,9),'r')
ylim([-0.25 1.5])
title ('PM')
subplot(2,3,5)
errorbar(TC_avg(:,5,1),TC_sem(:,5,1),'c')
hold on
errorbar(TC_avg(:,5,9),TC_sem(:,5,9),'r')
ylim([-0.25 1.5])
title ('AL')
subplot(2,3,6)
errorbar(TC_avg(:,6,1),TC_sem(:,6,1),'c')
hold on
errorbar(TC_avg(:,6,9),TC_sem(:,6,9),'r')
ylim([-0.25 1.5])
title ('LM')