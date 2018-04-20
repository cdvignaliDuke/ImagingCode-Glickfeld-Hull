mouse_mat = strvcat('i720', 'i738', 'i739');
area_mat = strvcat('V1','LM','PM','AL');
nmouse = size(mouse_mat,1);
narea = size(area_mat,1);
date_mat = cell(1,nmouse);
run_str_mat = cell(1,nmouse);
date_mat{1} = strvcat('170704','170628','170730','170802');
run_str_mat{1} = strvcat('runs-003-004','runs-002-003','runs-002-003','runs-002-003');
img_area_mat{1} = {'ML','PM','ML','V1'};
date_mat{2} = strvcat('170710','170712','170730','170804');
run_str_mat{2} = strvcat('runs-002-003','runs-002-003','runs-003-004','runs-002-003');
img_area_mat{2} = {'V1','LM','AL','PM'};
date_mat{3} = strvcat('170630','170711','170712','170731','170802','170803');
run_str_mat{3} = strvcat('runs-003-004','runs-002-003','runs-002-003','runs-002-003','runs-002-003','runs-002-003');
img_area_mat{3} = {'PM','ML','ML','ML','V1','ML'};

% date_mat.V1 = strvcat('170706','170710','170712');
% date_mat.LM = strvcat('170704','170712','170711');
% date_mat.PM = strvcat('170628', '170705', '170630');
% date_mat.AL = strvcat('170730', '170730', '170731');
% run_str_mat.V1 = strvcat('runs-003-004','runs-002-003','runs-002-003');
% run_str_mat.LM = strvcat('runs-003-004','runs-002-003','runs-002-003');
% run_str_mat.PM = strvcat('runs-002-003','runs-002-003','runs-003-004');
% run_str_mat.AL = strvcat('runs-002-003','runs-003-004','runs-002-003');

col_mat = strvcat('k','g','c','r');
frameRateHz = 30;
figure;

norm_all = struct;
[n,n2] = subplotn(nmouse);
for imouse = 1:nmouse
    mouse = mouse_mat(imouse,:);
    subplot(n,n2,imouse)
    area_str = [];
    for iarea = 1:size(area_mat,1)
        area_name = area_mat(iarea,:);
        ind = find(strcmp(img_area_mat{imouse},area_name));
        norm_temp = [];
        if length(ind)>0
            for i = 1:length(ind)
                date = eval(['date_mat{imouse}(' num2str(ind(i)) ',:)']);
                run_str = eval(['run_str_mat{imouse}(' num2str(ind(i)) ',:)']);
                load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PPresp.mat']))
                load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_decaySub.mat']))
                norm_off = bsxfun(@rdivide, resp_off, resp_off(:,end));
                norm_temp = [norm_temp; norm_off(good_ind,:)];
            end
        end
        ind2 = find(strcmp(img_area_mat{imouse}, 'ML'));
        if length(ind2)>0
            for i = 1:length(ind2)
                date = eval(['date_mat{imouse}(' num2str(ind2(i)) ',:)']);
                run_str = eval(['run_str_mat{imouse}(' num2str(ind2(i)) ',:)']);
                load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PPresp.mat']))
                load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_decaySub.mat']))
                load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_areaDivide.mat']))
                norm_off = bsxfun(@rdivide, resp_off, resp_off(:,end));
                area_ind = find(strcmp(cell_area,area_name));
                good_ind = intersect(good_ind,area_ind);
                norm_temp = [norm_temp; norm_off(good_ind,:)];
            end
        end
        if length(norm_temp)>1
            off_all = [off_all; 240];
            errorbar(off_all.*(1000/frameRateHz), squeeze(mean(norm_temp,1)),squeeze(std(norm_temp,[],1)./sqrt(size(norm_temp,1))), ['o' col_mat(iarea,:)]);
            hold on
            area_str = strvcat(area_str, [area_mat(iarea,:) ' n = ' num2str(size(norm_temp,1))]);
        else
            errorbar(off_all(1), 1,0,'ow')
            area_str = strvcat(area_str, [area_mat(iarea,:) ' n = 0']);
        end
        eval(['norm_all(imouse).' area_name '= norm_temp;']);
    end
    title(mouse)
    ylim([0 1.2])
    ylabel('Norm amplitude')
    xlabel('ISI (s)')
    legend(area_str,'location', 'southeast')
    axis square
end
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', ['ppAdaptation_eachAreaMouseSummary.pdf']),'-dpdf','-fillpage')

figure;
[n,n2] = subplotn(narea);
for iarea = 1:size(area_mat,1)
    area_name = area_mat(iarea,:);
    mouse_str = [];
    subplot(n,n2,iarea)
    for imouse = 1:nmouse
        mouse = mouse_mat(imouse,:);
        ind = find(strcmp(img_area_mat{imouse},area_name));
        norm_temp = [];
        if length(ind)>0
            for i = 1:length(ind)
                date = eval(['date_mat{imouse}(' num2str(ind(i)) ',:)']);
                run_str = eval(['run_str_mat{imouse}(' num2str(ind(i)) ',:)']);
                load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PPresp.mat']))
                load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_decaySub.mat']))
                norm_off = bsxfun(@rdivide, resp_off, resp_off(:,end));
                norm_temp = [norm_temp; norm_off(good_ind,:)];
            end
        end
        ind2 = find(strcmp(img_area_mat{imouse}, 'ML'));
        if length(ind2)>0
            for i = 1:length(ind2)
                date = eval(['date_mat{imouse}(' num2str(ind2(i)) ',:)']);
                run_str = eval(['run_str_mat{imouse}(' num2str(ind2(i)) ',:)']);
                load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PPresp.mat']))
                load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_decaySub.mat']))
                load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_areaDivide.mat']))
                norm_off = bsxfun(@rdivide, resp_off, resp_off(:,end));
                area_ind = find(strcmp(cell_area,area_name));
                good_ind = intersect(good_ind,area_ind);
                norm_temp = [norm_temp; norm_off(good_ind,:)];
            end
        end
        if length(norm_temp)>1
            off_all = [off_all; 240];
            errorbar(off_all.*(1000/frameRateHz), squeeze(mean(norm_temp,1)),squeeze(std(norm_temp,[],1)./sqrt(size(norm_temp,1))), ['o' col_mat(imouse,:)]);
            hold on
            mouse_str = strvcat(mouse_str, [mouse_mat(imouse,:) ' n = ' num2str(size(norm_temp,1))]);
        else
            errorbar(off_all(1), 1,0,'ow')
            mouse_str = strvcat(mouse_str, [mouse_mat(imouse,:) ' n = 0']);
        end
    end
    title(area_mat(iarea,:))
    ylim([0 1.2])
    legend(mouse_str,'location', 'southeast')
    axis square
    ylabel('Norm amplitude')
    xlabel('ISI (s)')
end
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', ['ppAdaptation_eachMouseAreaSummary.pdf']),'-dpdf','-fillpage')

figure;
subplot(1,2,1)
mouse_str = [];
for imouse = 1:nmouse
    temp = [];
    for iarea = 1:size(area_mat,1)
        temp = eval(['[temp; norm_all(imouse).' area_mat(iarea,:) ']']);
    end
    errorbar(off_all.*(1000/frameRateHz), squeeze(mean(temp,1)),squeeze(std(temp,[],1)./sqrt(size(temp,1))), ['o' col_mat(imouse,:)]);
    hold on
    mouse_str = strvcat(mouse_str, [mouse_mat(imouse,:) ' n = ' num2str(size(temp,1))]);
end
title('All areas')
ylim([0 1.2])
ylabel('Norm amplitude')
xlabel('ISI (s)')
axis square
legend(mouse_str,'location', 'southeast')

subplot(1,2,2)
area_str = [];
for iarea = 1:size(area_mat,1)
    temp = [];
    for imouse = 1:nmouse
        temp = eval(['[temp; norm_all(imouse).' area_mat(iarea,:) ']']);
    end
    errorbar(off_all.*(1000/frameRateHz), squeeze(mean(temp,1)),squeeze(std(temp,[],1)./sqrt(size(temp,1))), ['o' col_mat(iarea,:)]);
    hold on
    area_str = strvcat(area_str, [area_mat(iarea,:) ' n = ' num2str(size(temp,1))]);
end
title('All mice')
ylim([0 1.2])
ylabel('Norm amplitude')
xlabel('ISI (s)')
axis square
legend(area_str,'location', 'southeast')
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', ['ppAdaptation_allMouseAreaSummary.pdf']),'-dpdf','-fillpage')

%% orientation fits
mouse = mouse_mat(1,:);
date = date_mat{1}(1,:);
run_str = run_str_mat{1}(1,:);
load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
theta_smooth = 0:1:180;
k_hat = cell(nmouse,size(area_mat,1));
sse = cell(nmouse,size(area_mat,1));
sse_tot = cell(nmouse,size(area_mat,1));
R_square = cell(nmouse,size(area_mat,1));
y_fit = cell(nmouse,size(area_mat,1));
p = cell(nmouse,size(area_mat,1));
data_shift = cell(nmouse,size(area_mat,1));
ndir = length(dirs);
for iarea = 1:size(area_mat,1)
    fprintf([area_mat(iarea,:) '\r\n'])
    area_name = area_mat(iarea,:);
    mouse_str = [];
    for imouse = 1:nmouse
        fprintf([mouse_mat(imouse,:) '\r\n'])
        mouse = mouse_mat(imouse,:);
        ind = find(strcmp(img_area_mat{imouse},area_name));
        resp_temp = [];
        if length(ind>0)
            for i = 1:length(ind)
                date = eval(['date_mat{imouse}(' num2str(ind(i)) ',:)']);
                run_str = eval(['run_str_mat{imouse}(' num2str(ind(i)) ',:)']);
                load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PPresp.mat']))
                resp_temp = [resp_temp; resp_off_dir(good_ind,:,6)];
            end
        end
        ind2 = find(strcmp(img_area_mat{imouse},'ML'));
        if length(ind2>0)
            for i = 1:length(ind2)
                date = eval(['date_mat{imouse}(' num2str(ind2(i)) ',:)']);
                run_str = eval(['run_str_mat{imouse}(' num2str(ind2(i)) ',:)']);
                load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_areaDivide.mat']))
                load(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_PPresp.mat']))
                area_ind = find(strcmp(cell_area,area_name));
                good_ind = intersect(good_ind,area_ind);
                resp_temp = [resp_temp; resp_off_dir(good_ind,:,6)];
            end
        end
        if size(resp_temp,1)>0
            for iC = 1:size(resp_temp,1)
                data = resp_temp(iC,:);
                [val ind] = max(data,[],2);
                data_shift{imouse,iarea}(:,iC) = circshift(data,-ind+1)';
                [b_hat, k_hat{imouse,iarea}(1,iC), R_hat,u_hat,sse{imouse,iarea}(1,iC),R_square{imouse,iarea}(1,iC)] = miaovonmisesfit_ori(deg2rad(dirs),data);
                sse_tot{imouse,iarea}(1,iC) = -sse{imouse,iarea}(1,iC)./(R_square{imouse,iarea}(1,iC)- 1);
                fstat{imouse,iarea}(1,iC) = (sse{imouse,iarea}(1,iC)./(ndir-1))./(sse_tot{imouse,iarea}(1,iC)./(ndir-1));
                p{imouse,iarea}(1,iC) = 1-fcdf(fstat{imouse,iarea}(1,iC),ndir-1,ndir-1); 
                y_fit{imouse,iarea}(:,iC) = b_hat+R_hat.*exp(k_hat{imouse,iarea}(1,iC).*(cos(2.*(deg2rad(theta_smooth)-u_hat))-1));
            end
        end
    end
end

ntuned = zeros(1,size(area_mat,1));
ntot = zeros(1,size(area_mat,1));
osi = cell(1,size(area_mat,1));
for iarea = 1:size(area_mat,1)
    for imouse = 1:nmouse
        ind = find(p{imouse,iarea}>0.9);
        ntuned(1,iarea) = ntuned(1,iarea)+length(ind);
        ntot(1,iarea) = ntot(1,iarea)+size(p{imouse,iarea},2);
        [rpref prefind] = max(y_fit{imouse,iarea},[],1);
        orthind = prefind+90;
        orthind(find(orthind > 180)) = prefind(find(orthind > 180))-90;
        rorth = zeros(1,length(orthind));
        for iC = 1:length(orthind)
            rorth(:,iC) = y_fit{imouse,iarea}(orthind(iC),iC);
        end
        rorth(find(rorth<0)) = 0;
        osi{iarea} = [osi{iarea} (rpref(:,ind)-rorth(:,ind))./(rpref(:,ind)+rorth(:,ind))];
    end
end

figure;
cuts = [0.05 0.95 0.99 0.999];
for i = 1:4
    cut = cuts(i);
    subplot(2,2,i)
    area_str = [];
    for iarea = 1:size(area_mat,1)
        c = [];
        for imouse = 1:nmouse
            ind = find(p{imouse,iarea}>cut);
            c = [c data_shift{imouse,iarea}(:,ind)];
        end
        ind = find(c(1,:)>0);
        errorbar(dirs', mean(bsxfun(@rdivide,c(:,ind),c(1,ind)),2), std(bsxfun(@rdivide,c(:,ind),c(1,ind)),[],2)./sqrt(length(ind)),'-o')
        hold on
        area_str = strvcat(area_str, [area_mat(iarea,:) ' n = ' num2str(length(ind))]);
    end
    ylim([0 inf])
    xlim([0 180])
    title(num2str(cut))
    ylabel('Normalized amp')
    xlabel('Delta peak ori (deg)')
    legend(area_str,'Location', 'Northeast')
end
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', ['ppAdaptation_alignOriTuning.pdf']),'-dpdf','-fillpage')

figure;
for i = 1:4
    cut = cuts(i);
    subplot(2,2,i)
    ntuned = zeros(1,size(area_mat,1));
    ntot = zeros(1,size(area_mat,1));
    for iarea = 1:size(area_mat,1)
        for imouse = 1:nmouse
            ind = find(p{imouse,iarea}>cut);
            ntuned(1,iarea) = ntuned(1,iarea)+length(ind);
            ntot(1,iarea) = ntot(1,iarea)+size(p{imouse,iarea},2);
        end
    end
    [m, ci] = binofit(ntuned, ntot);
    errorbar(1:narea, m, m-ci(:,1)', ci(:,2)'-m, 'o')
    hold on
    ylim([0 1])
    xlim([0 narea+1])
    set(gca, 'XTick', 1:narea, 'XTickLabels', area_mat)
    axis square
    ylabel('Fraction tuned')
    title(num2str(cut))
end
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', ['ppAdaptation_fractionTuned.pdf']),'-dpdf','-fillpage')


figure;
for i = 1:4
    cut = cuts(i);
    subplot(2,2,i)
    osi = cell(1,size(area_mat,1));
    for iarea = 1:size(area_mat,1)
        for imouse = 1:nmouse
            ind = find(p{imouse,iarea}>cut);
            [rpref prefind] = max(y_fit{imouse,iarea},[],1);
            orthind = prefind+90;
            orthind(find(orthind > 180)) = prefind(find(orthind > 180))-90;
            rorth = zeros(1,length(orthind));
            for iC = 1:length(orthind)
                rorth(:,iC) = y_fit{imouse,iarea}(orthind(iC),iC);
            end
            rorth(find(rorth<0)) = 0;
            osi{iarea} = [osi{iarea} (rpref(:,ind)-rorth(:,ind))./(rpref(:,ind)+rorth(:,ind))];
        end
        errorbar(iarea, mean(osi{iarea},2), std(osi{iarea},[],2)./sqrt(size(osi{iarea},2)),'o');
        hold on
    end
    ylim([0 1])
    xlim([0 narea+1])
    set(gca, 'XTick', 1:narea, 'XTickLabels', area_mat)
    axis square
    ylabel('OSI')
    title(num2str(cut))
end
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', ['ppAdaptation_osi.pdf']),'-dpdf','-fillpage')

figure;
for i = 1:4
    cut = cuts(i);
    subplot(2,2,i)
    k = cell(1,size(area_mat,1));
    for iarea = 1:size(area_mat,1)
        for imouse = 1:nmouse
            ind = find(p{imouse,iarea}>cut);
            k{iarea} = [k{iarea} k_hat{imouse,iarea}(:,ind)];
        end
        errorbar(iarea, mean(k{iarea},2), std(k{iarea},[],2)./sqrt(size(k{iarea},2)),'o');
        hold on
    end
    ylim([0 3])
    xlim([0 narea+1])
    set(gca, 'XTick', 1:narea, 'XTickLabels', area_mat)
    axis square
    ylabel('K')
    title(num2str(cut))
end
print(fullfile('\\CRASH.dhe.duke.edu\data\home\lindsey\Analysis\2P', 'Adaptation', ['ppAdaptation_k.pdf']),'-dpdf','-fillpage')
