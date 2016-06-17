% date_mat = strvcat('150703', '150703','150704','150704','150720','150723', '151215'); %
% run1_mat = strvcat('_000_000', '_000_001', '_000_001', '_000_000','_000_001','_000_000', '_000_000'); %
% mouse_mat = strvcat('img24','img25','img24','img25','img28','img27', 'img32');   %
% subNum_mat = strvcat('924','925','924','925','928','927','932');   %
% time_mat = strvcat('1821', '2005','1845','1937','2015','1800','1423'); %
% nrun_mat = [1; 1; 1; 1; 1; 2; 1;]; %
% run2_mat = strvcat('XXXXXXXX', 'XXXXXXXX', 'XXXXXXXX', 'XXXXXXXX','XXXXXXXX','_000_001','XXXXXXXX'); %
% run_mat = cat(3,run1_mat,run2_mat);
% out_base = 'Z:\Analysis\2P Analysis\Lever\';
% base_dir =  'Z:\Data\2P_imaging\';

date_mat = strvcat('160517'); %
run1_mat = strvcat('_000_000'); %
mouse_mat = strvcat('img44');   %
subNum_mat = strvcat('944');   %
time_mat = strvcat('1601'); %
nrun_mat = [1]; %
run2_mat = strvcat('XXXXXXXX'); %
run_mat = cat(3,run1_mat,run2_mat);
out_base = 'Z:\Analysis\2P Analysis\Lever\';
base_dir =  'Z:\Data\2P_imaging\';


% %current datasets 11/9/15
% date_mat = strvcat('151016');%, '151018', '151029', '151102'); %
% run1_mat = strvcat('_000_000');%, '_000_000',  '_000_000', '_000_000'); %
% mouse_mat = strvcat('img30');%, 'img30', 'img30', 'img30');   %
% subNum_mat = strvcat('930');%, '930', '930', '930');   %
% time_mat = strvcat('1548');%,'1846','1127','1509'); %
% nrun_mat = [2]%; 2; 1; 2]; %
% run2_mat = strvcat('_000_001')%,'_000_001','XXXXXXXX','_000_001'); %
% run_mat = cat(3,run1_mat,run2_mat);
% out_base = 'Z:\Analysis\2P Analysis\Lever\';
% base_dir =  'Z:\Data\2P_imaging\';


%% all datasets
% date_mat = strvcat('150703', '150703','150704','150704','150720','150723', '151016', '151018', '151214', '151215', '160202','160203', '160319', '160320'); %
% run1_mat = strvcat('_000_000', '_000_001', '_000_001', '_000_000','_000_001','_000_000', '_000_000', '_000_000', '_000_000', '_000_000', '_000_000', '_000_002', '_000_001', '_000_001'); %
% mouse_mat = strvcat('img24','img25','img24','img25','img28','img27', 'img30', 'img30', 'img32', 'img32', 'img36', 'img36', 'img38', 'img38');   %
% subNum_mat = strvcat('924','925','924','925','928','927','930', '930', '932', '932', '936', '936', '938', '938');   %
% time_mat = strvcat('1821', '2005','1845','1937','2015','1800','1548','1846', '1423', '1655', '1531','1535', '1723', '1347'); %
% nrun_mat = [1; 1; 1; 1; 1; 2; 2; 2; 1; 1; 1; 1; 1; 1;]; %
% run2_mat = strvcat('XXXXXXXX', 'XXXXXXXX', 'XXXXXXXX', 'XXXXXXXX','XXXXXXXX','_000_001','_000_001','_000_001', 'XXXXXXXX', 'XXXXXXXX', 'XXXXXXXX', 'XXXXXXXX', 'XXXXXXXX', 'XXXXXXXX'); %
% run_mat = cat(3,run1_mat,run2_mat);
% out_base = 'Z:\Analysis\2P Analysis\Lever\';
% base_dir =  'Z:\Data\2P_imaging\';
% % 

% date_mat = strvcat('151215', '160202','160203', '160319', '160320', '160324');%, 
% run1_mat = strvcat('_000_000', '_000_000', '_000_002', '_000_001', '_000_000', '_000_001');%,
% mouse_mat = strvcat('img32', 'img36', 'img36', 'img38', 'img38', 'img41');%,
% subNum_mat = strvcat('932', '936', '936', '938', '938', '941');%
% time_mat = strvcat('1655', '1531','1535', '1723', '1347', '1450');
% nrun_mat = [1; 1; 1; 1; 1; 1];
% if nrun_mat == 1;
%     run_mat =run1_mat; 
% else
%     run2_mat = strvcat('_000_001')
%     run_mat = cat(3,run1_mat,run2_mat);
% end
% out_base = 'Z:\Analysis\2P Analysis\Lever\';
% base_dir =  'Z:\Data\2P_imaging\';

% 
% 
% date_mat = strvcat('151215');%, 
% run1_mat = strvcat('_000_000');%,
% mouse_mat = strvcat('img32');%,
% subNum_mat = strvcat('9432');%
% time_mat = strvcat('1423');
% nrun_mat = [1;];
% if nrun_mat == 1;
%     run_mat =run1_mat; 
% else
%     run2_mat = strvcat('_000_001')
%     run_mat = cat(3,run1_mat,run2_mat);
% end
% out_base = 'Z:\Analysis\2P Analysis\Lever\';
% base_dir =  'Z:\Data\2P_imaging\';
