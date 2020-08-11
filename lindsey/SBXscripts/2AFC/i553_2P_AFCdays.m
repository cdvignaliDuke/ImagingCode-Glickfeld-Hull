mouse = 'i553';
date = strvcat('170216','170217','170220','170221','170222','170223','170224','170226','170227','170228','170301','170302','170303');
run = [{'001'}, {'001'}, {'001'}, {'001'}, {'001'}, {'001'}, {'001'}, {'001'}, {'001'}, {'001'}, {'001'}, {'003'}, {strvcat('001', '002')}];
nrun = [1,1,1,1,1,1,1,1,1,1,1,1,2];
data_avg_mat = [];
dfof_max_mat = [];
for i = 1:size(date,1)
    run_str = catRunName(run{i}, nrun(:,i));
    if exist(fullfile(['\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P'], [date(i,:) '_' mouse], [date(i,:) '_' mouse '_' run_str], [date(i,:) '_' mouse '_' run_str '_reg_shifts.mat']))
        load(fullfile(['\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P'], [date(i,:) '_' mouse], [date(i,:) '_' mouse '_' run_str], [date(i,:) '_' mouse '_' run_str '_reg_shifts.mat']));
        data_avg_mat = cat(3,data_avg_mat,data_avg);
    else
        data_avg_mat = cat(3,data_avg_mat,NaN(size(data_avg)));
    end
    if exist(fullfile(['\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P'], [date(i,:) '_' mouse], [date(i,:) '_' mouse '_' run_str], [date(i,:) '_' mouse '_' run_str '_mask_cell.mat']))
        load(fullfile(['\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P'], [date(i,:) '_' mouse], [date(i,:) '_' mouse '_' run_str], [date(i,:) '_' mouse '_' run_str '_mask_cell.mat']))
        dfof_max_mat = cat(3,dfof_max_mat,data_dfof_max);
    else
        dfof_max_mat = cat(3,dfof_max_mat,NaN(size(data_dfof_max)));
    end
end

figure;
[n n2] =subplotn(size(date,1));
for i = 1:size(date,1)
    subplot(n,n2,i)
    imagesc(data_avg_mat(:,:,i))
    title(date(i,:))
end

figure;
[n n2] =subplotn(size(date,1));
for i = 1:size(date,1)
    subplot(n,n2,i)
    imagesc(dfof_max_mat(:,:,i))
    title(date(i,:))
end