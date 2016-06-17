%plot corrCoef between df/f and licking
colors = {'b', 'b', 'r', 'r', 'y', 'c', 'c', 'k', 'b', 'b', 'r', 'm', 'r'};
days = {'160319_img41', '160320_img41', '160314_img38', '160315_img38', '160129_img36', '160129_img35', '160131_img35', '151211_img32', '150717_img28', '150716_img28', '150718_img27', '151022_img29', '150719_img27'};
DATA_DIR = 'Z:\Analysis\LeverAnalysis\CorrCoefSummary\';
%DATA_DIR = 'Z:\Analysis\LeverAnalysis\LeverSummaryNoShift\';
summary_coef = {}; 
for kk = 1:length(days)
    curr_file_coef = strcat(DATA_DIR, days{kk}, '_corrCoefLick');
    summary_coef{kk} = load(curr_file_coef);
end 
coefMat = [];
for i = 1:length(summary_coef)
coefMat = [coefMat, mean(summary_coef{i}.coefMat2)];
end
figure;
for i = 1:length(coefMat)
    if i <5
        plot(1,coefMat(i), ['x' colors{i}]); hold on;
    else
        plot(1,coefMat(i), ['o' colors{i}], 'MarkerFaceColor', colors{i}); hold on;
    end
end
hline(0,'--k')
vline(0,'--k')
ylim([-1 1])
xlim([.5 2])
legend(days{1:length(days)})