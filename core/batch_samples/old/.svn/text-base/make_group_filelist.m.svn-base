
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filelist='exp_list_scmr8_3.xls';
s=readXLSfilelist(filelist);

groups=[s.group];
Ngroups=max(groups);

for i=1:Ngroups
    run_index=find(groups==i);
    if isempty(run_index)
        continue;
    end
    
    groupname=['group',num2str(i)];

    for j=1:length(run_index)
        s(run_index(j)).filename=[s(run_index(j)).filename,'_',groupname];
    end
end

s(find(isnan(groups)))=[];

