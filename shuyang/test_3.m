%% topology - optimal resolution - test

resolution_ratios  = [5,250,500,750,1000,1250,1500,1750,2000,20000]; %divide the total process into how many parts
speed_topo_ave = zeros(1,length(resolution_ratios));
for r = 1: length(resolution_ratios)
    speed_topo_mat= zeros(length(speed_run_cell),resolution_ratios(r));
    for i=1:1:length(speed_run_cell)
        speed_temp = speed_run_cell{i};
        mapping = linspace(1,length(speed_temp)+1,resolution_ratios(r));
        for j = 1:1:length(speed_temp)
            speed_topo_mat(i,(mapping>=j))=speed_temp(j);
        end
    end
    speed_topo_ave(r) = std(mean(speed_topo_mat));
end
figure;plot(speed_topo_ave);


