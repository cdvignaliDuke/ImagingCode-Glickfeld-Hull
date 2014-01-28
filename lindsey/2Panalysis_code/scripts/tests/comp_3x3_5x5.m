P = 1;
matrix = 'SF3xTF3';
image = 'axon';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
base = 'G:\users\lindsey\analysisLG\active mice';
post_win = [5 14];

fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_POST' num2str(post_win) 'fit_all_summary.mat']);
load(fn);

Fit_all_3x3 = Fit_all_invlog;
Speed_3x3 = Speed;
H_ttest_all_3x3 = H_ttest_all;
mice_3x3 = mice;
mice_list_3x3 = strvcat(mice_3x3);
nexp_per_area_3x3= nexp_per_area;

P = 1;
matrix = 'SF5xTF5';
image = 'axon';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
base = 'G:\users\lindsey\analysisLG\active mice';
post_win = [5 14];

fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_POST' num2str(post_win) '_fit_all_summary.mat']);
load(fn);

Fit_all_5x5 = Fit_all_invlog;
Speed_5x5 = Speed;
H_ttest_all_5x5 = H_ttest_all;
mice_5x5 = mice;
mice_list_5x5 = strvcat(mice_5x5);
nexp_per_area_5x5= nexp_per_area;

P = 1;
matrix = 'SF5xTF5';
image = 'axon';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
base = 'G:\users\lindsey\analysisLG\active mice';
post_win = [5 14];

fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_POST' num2str(post_win) '_subset_fit_all_summary.mat']);
load(fn);

Fit_all_5x5_subset = Fit_all_invlog;
Speed_5x5_subset = Speed;
H_ttest_all_5x5_subset = H_ttest_all;
mice_5x5_subset = mice;
mice_list_5x5_subset = strvcat(mice_5x5);
nexp_per_area_5x5_subset= nexp_per_area;

ind = zeros(1,length(mice_list_5x5));
for mouse5 = 1:length(mice_list_5x5);
    for mouse3 = 1:length(mice_list_3x3);
        if mice_list_5x5(mouse5,:) == mice_list_3x3(mouse3,:);
            ind(:,mouse5) = mouse3;
        end
    end
end

figure;
subplot(3,3,1);
for mouse = 1:length(mice_list_5x5);
    for iArea = 1:8;
        if H_ttest_all_5x5(iArea,1,mouse)==1;
        	if H_ttest_all_3x3(iArea,1,ind(:,mouse))==1;
                plot(Fit_all_3x3(iArea,4,1,ind(:,mouse)),Fit_all_5x5(iArea,4,1,mouse),['*' col(iArea)])
                hold on
            end
        end
    end
end
x = [0:.01:.35];
y= x;
plot(x,y,'-k')
xlabel('3x3')
ylabel('5x5')
title('SF');

subplot(3,3,2);
for mouse = 1:length(mice_list_5x5);
    for iArea = 1:8;
        if H_ttest_all_5x5(iArea,1,mouse)==1;
        	if H_ttest_all_3x3(iArea,1,ind(:,mouse))==1;
                plot(Fit_all_3x3(iArea,5,1,ind(:,mouse)),Fit_all_5x5(iArea,5,1,mouse), ['*' col(iArea)])
                hold on
            end
        end
    end
end
x = [0:1:20];
y= x;
plot(x,y,'-k')
xlabel('3x3')
ylabel('5x5')
title('TF');

subplot(3,3,3);
for mouse = 1:length(mice_list_5x5);
    for iArea = 1:8;
        if H_ttest_all_5x5(iArea,1,mouse)==1;
        	if H_ttest_all_3x3(iArea,1,ind(:,mouse))==1;
                plot(Speed_3x3(iArea,1,ind(:,mouse)),Speed_5x5(iArea,1,mouse), ['*' col(iArea)])
                hold on
            end
        end
    end
end
x = [0:1:800];
y= x;
plot(x,y,'-k')
xlim([0 800])
ylim([0 800])
xlabel('3x3')
ylabel('5x5')
title('Speed');

subplot(3,3,4);
for mouse = 1:length(mice_list_5x5);
    for iArea = 1:8;
        if H_ttest_all_5x5_subset(iArea,1,mouse)==1;
        	if H_ttest_all_3x3(iArea,1,ind(:,mouse))==1;
                plot(Fit_all_3x3(iArea,4,1,ind(:,mouse)),Fit_all_5x5_subset(iArea,4,1,mouse),['*' col(iArea)])
                hold on
            end
        end
    end
end
x = [0:.01:.35];
y= x;
plot(x,y,'-k')
xlabel('3x3')
ylabel('5x5 subset')
title('SF');

subplot(3,3,5);
for mouse = 1:length(mice_list_5x5);
    for iArea = 1:8;
        if H_ttest_all_5x5_subset(iArea,1,mouse)==1;
        	if H_ttest_all_3x3(iArea,1,ind(:,mouse))==1;
                plot(Fit_all_3x3(iArea,5,1,ind(:,mouse)),Fit_all_5x5_subset(iArea,5,1,mouse), ['*' col(iArea)])
                hold on
            end
        end
    end
end
x = [0:1:20];
y= x;
plot(x,y,'-k')
xlabel('3x3')
ylabel('5x5 subset')
title('TF');

subplot(3,3,6);
for mouse = 1:length(mice_list_5x5);
    for iArea = 1:8;
        if H_ttest_all_5x5_subset(iArea,1,mouse)==1;
        	if H_ttest_all_3x3(iArea,1,ind(:,mouse))==1;
                plot(Speed_3x3(iArea,1,ind(:,mouse)),Speed_5x5_subset(iArea,1,mouse), ['*' col(iArea)])
                hold on
            end
        end
    end
end
x = [0:1:800];
y= x;
plot(x,y,'-k')
xlim([0 800])
ylim([0 800])
xlabel('3x3')
ylabel('5x5 subset')
title('Speed');

subplot(3,3,7);
for mouse = 1:length(mice_list_5x5);
    for iArea = 1:8;
        if H_ttest_all_5x5(iArea,1,mouse)==1;
            plot(Fit_all_5x5(iArea,4,1,mouse),Fit_all_5x5_subset(iArea,4,1,mouse),['*' col(iArea)])
            hold on
        end
    end
end
x = [0:.01:.35];
y= x;
plot(x,y,'-k')
xlabel('5x5')
ylabel('5x5 subset')
title('SF');

subplot(3,3,8);
for mouse = 1:length(mice_list_5x5);
    for iArea = 1:8;
        if H_ttest_all_5x5(iArea,1,mouse)==1;
            plot(Fit_all_5x5(iArea,5,1,mouse),Fit_all_5x5_subset(iArea,5,1,mouse), ['*' col(iArea)])
            hold on
        end
    end
end
x = [0:1:20];
y= x;
plot(x,y,'-k')
xlabel('5x5')
ylabel('5x5 subset')
title('TF');

subplot(3,3,9);
for mouse = 1:length(mice_list_5x5);
    for iArea = 1:8;
        if H_ttest_all_5x5(iArea,1,mouse)==1;
            plot(Speed_5x5(iArea,1,mouse),Speed_5x5_subset(iArea,1,mouse), ['*' col(iArea)])
            hold on
        end
    end
end
x = [0:1:800];
y= x;
plot(x,y,'-k')
xlim([0 800])
ylim([0 800])
xlabel('5x5')
ylabel('5x5 subset')
title('Speed');

fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_POST' num2str(post_win) '5x5_vs_3x3_vs_5x5subset.ps']);
print(gcf, '-depsc', fn_out);