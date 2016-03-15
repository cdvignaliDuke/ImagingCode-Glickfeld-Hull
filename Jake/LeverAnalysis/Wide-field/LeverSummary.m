%% SUMMARY SCATTER PLOT
colors = {'r', 'r', 'r', 'r', 'r', 'b', 'b', 'b', 'b', 'm', 'm', 'm', 'g', 'g', 'g'};
days = {'150514_img24', '150518_img24', '150519_img24', '150521_img24' '150706_img24', '150517_img25', '150518_img25', '150514_img25', '150515_img25', '150716_img27', '150718_img27', '150719_img27', '150716_img28', '150717_img28', '150719_img28'};
DATA_DIR = 'C:\Users\jake\TempData\summaryFolder\';
summary_succ = {}; 
summary_fail = {}; 
for kk = 1:length(days)
    curr_file_succ = strcat(DATA_DIR, days{kk}, '_success');
    summary_succ{kk} = load(curr_file_succ);
    
    curr_file_fail = strcat(DATA_DIR, days{kk}, '_fail');
    
    temp2 = load(curr_file_fail);
    summary_fail{kk} = temp2;
    
end 

fail_avg = [];
succ_avg = []; 
fail_std = [];
succ_std = []; 
fail_sm = []; 
succ_sm = []; 

for kk = 1:length(days)
succ_avg(1,kk) = mean(mean(summary_succ{kk}.success_roi(:,7:9),1),2);
fail_avg(1,kk) = mean(mean(summary_fail{kk}.fail_roi(:,7:9),1),2);
end
for kk = 1:length(days)
succ_sm(1,kk) = std(mean(summary_succ{kk}.success_roi(:,7:9),2),[],1)./sqrt(size(summary_succ{kk}.success_roi,1));
fail_sm(1,kk) = std(mean(summary_fail{kk}.fail_roi(:,7:9),2),[],1)./sqrt(size(summary_fail{kk}.fail_roi,1));
end

figure;
for i = 1:length(succ_avg); 
    errorbarxy(succ_avg(i), fail_avg(i), succ_sm(i), fail_sm(i), [],[],'w', colors{i}); hold on; 
    plot(succ_avg(i), fail_avg(i), ['o' colors{i}]); hold on;
end
ylim([-.03 .15])
xlim([-.03 .15])
x = -.1:.1:1;
y = x;
hold on; plot(x,y,'k')
hline(0,'--k')
vline(0,'--k')

errorbarxy(fail_avg(i), succ_avg(i), fail_sm(i), succ_sm(i), [],[],'w', colors{i}); hold on; 
    plot(fail_avg(i), succ_avg(i), ['o' colors{i}]); hold on;
    

%% spatial correlation scatter plot 
img24 = [0.60, 0.81, 0.47];
img25 = [0.93, 0.96, 0.80, 0.92, 0.79, 0.80];
img27 = [0.76, 0.65, 0.64];
%img28 = 

valuesY = [img24, img25, img27];
categoriesX = [1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3];

corr = scatter( categoriesX, valuesY, 20, 'b', 'filled')
xlim([0 4]);
ylim([0 1]);
corr.XTick = [0 1 2 3 4];

    
    
    