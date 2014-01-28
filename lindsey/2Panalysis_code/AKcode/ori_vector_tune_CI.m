function [LB68 LB95 UB68 UB95 CI68 CI95 boot_vector_tune]=ori_vector_tune_CI(ori_table)

% si=size(dir_table);
% ori_table=reshape(dir_table,si(1),si(2)./2,2);
% ori_table=squeeze(nanmean(ori_table,3));

si=size(ori_table);
ori_resp=squeeze(nanmean(ori_table));

[angle_avg,mag_avg,tune_avg]...
      = vector_average(ori_resp);
  
for i=1:1000
    for j=1:si(2)
        ori_table_rand(:,j)=randsample(ori_table(:,j),si(1),'true');
    end
    ori_resp_rand = squeeze(nanmean(ori_table_rand));
    [boot_vector_tune(i) mag] = vector_tune_fixed(ori_resp_rand,angle_avg);
end


LB68=prctile(boot_vector_tune,15.9);
UB68=prctile(boot_vector_tune,84.1);
LB95=prctile(boot_vector_tune,2.5);
UB95=prctile(boot_vector_tune,97.5);
CI68=UB68-LB68;
CI95=UB95-LB95;