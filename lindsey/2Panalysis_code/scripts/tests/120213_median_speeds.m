for iArea = 1:3
    nexp = all_fits(iArea).nexp;
    for iexp = 1:nexp
        if all_fits(iArea).expt(iexp).n(2)>0
            n = all_fits(iArea).expt(iexp).n(1);
            speed = [];
            for iCell = 1:n
                if all_fits(iArea).expt(iexp).bouton(iCell).goodfit ==1
                    speed = [speed all_fits(iArea).expt(iexp).bouton(iCell).speed];
                end
            end
            all_fits(iArea).expt(iexp).median_speed = median(speed);
            all_fits(iArea).expt(iexp).mean_speed = mean(speed);
        end
    end
end


area_order = [2; 3; 1];
avg_median_speed = zeros(3,2);
for iArea = 1:3
    nexp = all_fits(area_order(iArea)).nexp;
    median_speed = [];
    for iexp = 1:nexp
        if all_fits(area_order(iArea)).expt(iexp).n(2)>50
            median_speed = [median_speed all_fits(area_order(iArea)).expt(iexp).median_speed];
            hold on
        end
    end
    avg_median_speed(iArea,1) = mean(median_speed);
    avg_median_speed(iArea,2) = std(median_speed,[],2)./sqrt(length(median_speed));
end


mouse_list = {'Y13' 'X32' 'DR7' 'DR9' 'AC39' 'AC42' 'AC44' 'AC45' 'Y18' 'Y26' 'M13' 'M14' 'M22' 'M31'};
inj = 'V1';
matrix = 'SF5xTF5';
P = 2;
cutoff = [0 10 25 50];
figure;
for icut = 1:4
    subplot(2,2,icut)
for iMouse = 1:length(mouse_list);
    mouse = mouse_list{iMouse};
    medians = [NaN; NaN; NaN];
    for iArea = 1:3
        image = areas(area_order(iArea),:);
        sum_base = 'G:\users\lindsey\analysisLG\experiments';
        list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
        load(list_fn);
        median_speed = [];
        nexp = all_fits(area_order(iArea)).nexp;
        for iexp = 1:nexp
            exp_mouse = exp_list.mouse_mat{iexp};
            if length(exp_mouse)==length(mouse)
                if exp_mouse == mouse;
                    if all_fits(area_order(iArea)).expt(iexp).n(2)>cutoff(icut)
                        median_speed = [median_speed all_fits(area_order(iArea)).expt(iexp).median_speed];
                    end
                end
            end
        end
        if length(median_speed)>0
            if length(median_speed)>1
                median_speed = mean(median_speed);
            end
            medians(iArea,:) = median_speed;
        end
    end
    semilogy(1:1:3, medians, '-ob')
    hold on
end
    title(['at least ' num2str(cutoff(icut)+1) ' boutons'])
    xlim([0 4])
end

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_median_speed_allmice.pdf']);
    print(gcf, '-dpdf', fn_out);

figure;
area_order = [1 2; 2 3; 1 3];
for iMouse = 1:length(mouse_list);
    mouse = mouse_list{iMouse};
    for ipair = 1:3
        medians = [NaN NaN];
        for iArea = 1:2
            image = areas(area_order(ipair, iArea),:);
            sum_base = 'G:\users\lindsey\analysisLG\experiments';
            list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
            load(list_fn);
            median_speed = [];
            nexp = all_fits(area_order(ipair,iArea)).nexp;
            for iexp = 1:nexp
                exp_mouse = exp_list.mouse_mat{iexp};
                if length(exp_mouse)==length(mouse)
                    if exp_mouse == mouse;
                        if all_fits(area_order(ipair,iArea)).expt(iexp).n(2)>25
                            median_speed = [median_speed all_fits(area_order(ipair,iArea)).expt(iexp).median_speed];
                        end
                    end
                end
            end
            if length(median_speed)>0
                if length(median_speed)>1
                    median_speed = mean(median_speed);
                end
                medians(:,iArea) = median_speed;
            end
        end
        subplot(3,3,ipair+3)
        semilogy(1:2, medians, '-ob')
        hold on
        title([areas(area_order(ipair,1),:) ' vs ' areas(area_order(ipair,2),:)])
        xlim([0 3])
        ylim([1 1000])
    end    
end

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_median_speed_allmice_areapairs.pdf']);
    print(gcf, '-dpdf', fn_out);
