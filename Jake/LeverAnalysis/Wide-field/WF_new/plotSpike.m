

nDraws = nIC;
sampleN = randperm(nIC, nDraws);
sampleN = 1:nIC;

for ic = 1:nDraws%57 used for 2p example  cell 11 for img38
%     subplot(2,2,ic)
    tt = (1:size(data_tc_spont(:,sampleN(ic)),1)) * 100;
    figure;
    plot(tt, data_tc_spont(:,sampleN(ic)) - 20, 'k');
%     plot( data_tc_spont(:,sampleN(ic)) + 2, 'k');
    hold on;
    plot(tt(1:end-1), events_diff_x1(:,sampleN(ic)), 'b');
%     plot( events_diff_x1(:,sampleN(ic)), 'b');
    scale = max(events_diff_x1(:,sampleN(ic)));
    
    plot( events_ind{sampleN(ic)}*100, scale*0.6, 'k.');
%     plot( events_ind{sampleN(ic)}, scale*0.6, 'k.');
        
%     for k = 1:nIC
%     %     events_ind{ic} = find(events_diff_x(:,ic)>thresh(ic));
%     [~, events_ind{k}, ~] = CellsortFindspikes(events_diff_x1(:,k), thresh, dt, deconvtau, normalization);
%     events_ind{k}(find(diff(events_ind{k})==1)+1)=[];
% %     events_rate(1,ic) = length(events_ind{ic})./((double(ifi)./1000)*size(data_tc_spont,1));
%     end
%     plot(events_diff_x2(:,sampleN(ic)) - 10, 'm');
%     scale = max(events_diff_x1(:,sampleN(ic)));
%     plot(events_ind{sampleN(ic)} , scale*0.6 - 10, 'r*');
%     xlim([3000,3500]);
end
% thresh = 2.5;
title('Raw Trace--Black  1st Derivative--Blue');
supertitle([date, ' Threshold = ', num2str(thresh), '  1st Derivative--Black', '   2nd Derivative--Purple']);
saveas(fig,[data_dir, date, 'std', num2str(thresh), 'spike_threshold_2ndDev.fig']);
