r1 = load('Z:\home\jake\Analysis\Cue_reward_pairing_analysis\2P\170420_001_img92\thresholding\_cell_TCs_half.mat');
r2 = load('Z:\home\jake\Analysis\Cue_reward_pairing_analysis\2P\170427_000_img92\thresholding\_cell_TCs.mat');
avg_all_nolick = [r1.avg_NR_nolick r1.avg_OR_nolick r1.avg_UR_nolick];
avg_all_nolick2 = [r2.avg_NR_nolick r2.avg_OR_nolick r2.avg_UR_nolick];
ymax = max(max(avg_all_nolick,[],2),[],1);
ymin = min(min(avg_all_nolick,[],2),[],1);

figure;
i = 1;
for ic = [10 26 99];
    subplot(3,2,(i-1)*2+1);
    plot(tt,r1.avg_NR_nolick(ic,:),'k');hold on;
    
    plot(tt,r1.avg_OR_nolick(ic,:),'r');
    ylim([ymin*1.1 ymax*1.1])
    xlim([-1 2])
    i = i+1;
end

i = 1;
for ic = [26 33 82];
    subplot(3,2,(i-1)*2+2);
    plot(tt,r2.avg_NR_nolick(ic,:),'k');hold on;
    
    plot(tt,r2.avg_OR_nolick(ic,:),'r');
    ylim([ymin*1.1 ymax*1.1])
    xlim([-1 2])
    i = i+1;
end
% if ~isempty(UR_movie_nolick)
%     plot(tt,avg_UR_nolick(ic,:),'g');
% end

suptitle([date ' ' mouse ' Average cue resp: Normal- black (n = ' num2str(size(NR_movie,1)) ' trials); Omit- red (n = ' num2str(size(OR_movie,1)) ' trials); Unexpect- green (n = ' num2str(size(UR_movie,1)) ' trials)'])
orient landscape