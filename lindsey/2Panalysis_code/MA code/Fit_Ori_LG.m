%Fitting ori x SF: 

%PLAN:
%0) plot raw data, ori vs SF
%1) average 2 dirns, then find best ori (avg'ed over SF).
%2) extract SF params using DOG, using best ori
%3) separately, find best SF (at pref ori), then fit ori tuning curve at that
%best SF
%4) repeat steps 2-3 100 times after sampling with replacement, get bootstrap
%estimates for SF, SF_width, ori and ori_width

s = struct;
s.data = data;
ind0 = find(data<0);
data1 = data;
data1(ind0) = 0;
s.data1 = data1;
s.data_std = data_std;
s.dF = dF_mat(iCell,:)';
s.ypos = i(iCell, :);
s.xpos = j(iCell, :);
%first, collapse across directions: 
data2 = (data(1:dirs/2,:) + data(1+(dirs/2):dirs,:))./2;
ind0 = find(data2<0);
data2(ind0) = 0;
s.data2 = data2;
s.y = ORI_vec0;
s.sfsf_tftf_grid = [sfsf_grid; tftf_grid];
%step 1: 
%find approx best ori: 
if nSFTF > 1;
    data3 = mean(data2,2);
    [max_resp, ind_bestori] = max(data3,[],1);
    s.ind_bestori = ind_bestori;
    s.data3 = data3;

    %step2:
    %now choose the peak SF and fit ori tuning curve (and then estimate
    %direction tuning..

    [max_resp, ind_bestSFTF] = max(data2(ind_bestori,:));
    s.ind_bestSFTF = ind_bestSFTF;
    s.best_SF = sfsf_grid(1,ind_bestSFTF);
    s.best_TF = tftf_grid(1,ind_bestSFTF);
elseif nSFTF == 1
    [max_resp, ind_bestori] = max(data2,[],1);
    s.ind_bestori = ind_bestori;
    s.data3 = data2;
    s.best_SF = SFs;
    s.best_TF = TFs;
    ind_bestSFTF = 1;
    s.ind_bestSFTF = 1;
end

%now  compute classical ori estimates

for iSFTF = 1:nSFTF
    OriStatKO.ori_ratio_change = data2(:,iSFTF);
    OriStatKO.dir_ratio_change = data1(:,iSFTF);
    % classical analysis
    [OriStatKO.best_dir,...
        OriStatKO.null_dir,...
        OriStatKO.R_best_dir,...
        OriStatKO.R_null_dir,...
        OriStatKO.R_min_dir,...
        OriStatKO.DI]...
        = dir_indexKO(OriStatKO.dir_ratio_change);
    [OriStatKO.best_ori,...
        OriStatKO.null_ori,...
        OriStatKO.R_best_ori,...
        OriStatKO.R_null_ori,...
        OriStatKO.R_min_ori,...
        OriStatKO.OI]...
        = dir_indexKO(OriStatKO.ori_ratio_change);

    % vector average analysis

    [OriStatKO.dir_vector_angle,...
        OriStatKO.dir_vector_mag,...
        OriStatKO.dir_vector_tune]...
        = vector_average(OriStatKO.dir_ratio_change);

    [OriStatKO.ori_vector_angle,...
        OriStatKO.ori_vector_mag,...
        OriStatKO.ori_vector_tune]...
        = vector_average(OriStatKO.ori_ratio_change);

    OriStatKO.ori_vector_angle = OriStatKO.ori_vector_angle /2;

    % alternate tuning measure

    OI2 = (OriStatKO.R_best_ori - OriStatKO.R_null_ori)./(OriStatKO.R_best_ori + OriStatKO.R_null_ori);
    DI2 = (OriStatKO.R_best_dir - OriStatKO.R_null_dir)./(OriStatKO.R_best_dir + OriStatKO.R_null_dir);


    s.OriStatKO(iSFTF,:) = OriStatKO;
    s.x(iSFTF,:) = [OriStatKO.OI OriStatKO.ori_vector_tune OI2 OriStatKO.best_ori OriStatKO.ori_vector_angle...
        OriStatKO.DI OriStatKO.dir_vector_tune DI2 OriStatKO.best_dir OriStatKO.dir_vector_angle];
end

if SAVEALLDATA == 0
    s = rmfield(s,'ind_bestori');
    s = rmfield(s,'data');
    s = rmfield(s,'data_std');
    s = rmfield(s,'data1');
    s = rmfield(s,'data2');
    s = rmfield(s,'data3');
    s = rmfield(s,'best_TF');
    s = rmfield(s,'best_SF');
    s = rmfield(s,'OriStatKO');
    s = rmfield(s,'dF');
    s = rmfield(s,'xpos');
    s = rmfield(s,'ypos');
    s = rmfield(s,'sfsf_tftf_grid');
end
