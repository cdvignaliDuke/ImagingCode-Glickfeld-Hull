%plot_session_PSTH
% plot the mean PSTH for a given session
% also plot the mean PSTH for each neuron
function plot_session_PSTH(normalR_hist_D1, omitR_hist_D1, plot_ops, pdf_dir)
x_axis = plot_ops.x_axis;
%plot mean across cells and trials
fig1 = figure('rend', 'painters', 'pos', [50 150 750 550]); 
errorbar(x_axis, mean(normalR_hist_D1,2), std(normalR_hist_D1,[],2)/sqrt(size(normalR_hist_D1,2)), 'k'); hold on;
errorbar(x_axis, mean(omitR_hist_D1,2), std(omitR_hist_D1,[],2)/sqrt(size(omitR_hist_D1,2)), 'r');
plot(x_axis, mean(normalR_hist_D1,2), 'k');
plot(x_axis, mean(omitR_hist_D1,2), 'r');
xlim([-500 2000]);
ylim([0 3]);
vline(plot_ops.cue, 'c');
vline(plot_ops.reward, 'b');
ylabel('Firing rate (Hz)');
xlabel('time from cue (ms)');
title([plot_ops.mouse, ' ', plot_ops.day, ' mean PSTH. n=', num2str(size(normalR_hist_D1,2))]);
print([pdf_dir, plot_ops.mouse, ' ', plot_ops.day, ' mean PSTH.pdf'], '-dpdf');

%plot mean across trials for each cell
cell_n = size(normalR_hist_D1,2);
fig_n = ceil(cell_n/42);
y_lim = [0 max(max([normalR_hist_D1; omitR_hist_D1]))];
this_cell = 1;
for this_fig = 1:length(fig_n)
    figure('rend', 'painters', 'pos', [50 50 750 750]);
    for this_subplot = 1:42
        subplot(7,6,this_subplot);
        plot(x_axis, omitR_hist_D1(:,this_cell), 'r');  hold on;
        plot(x_axis, normalR_hist_D1(:,this_cell), 'k');
        this_cell = this_cell +1;
        ylim([0 3]);
        xlim([-500 1500]);
        vline(plot_ops.cue, 'c');
        vline(plot_ops.reward, 'b');
        if this_cell == cell_n
            break
        end
    end
    suptitle([plot_ops.mouse, ' ', plot_ops.day, ' all cells ', num2str(this_fig)]);
    print([pdf_dir, plot_ops.mouse, ' ', plot_ops.day, ' all cells ', num2str(this_fig), '.pdf'], '-dpdf');
end
return 