pairs=[[1,2];
        [19,20];
        [25,26];
        [35,36];
        [37,38];
        [39,40];
        [41,42];
        [43,44];
        [47,53];
        [48,53];];

filelist='exp_list_scmr6_2.xls';
data_dir='\\zmey\storlab\data\Soumya\scmr6\';
out_dir='E:\users\kenichi\temp';
    
for i=1:size(pairs,1)    
    analyze_2stacks_stats_oristats (filelist, pairs(i,:), data_dir, out_dir);
end
    
    
        
        
    