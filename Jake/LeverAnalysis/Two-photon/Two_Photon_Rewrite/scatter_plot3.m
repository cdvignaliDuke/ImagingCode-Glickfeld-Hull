function scatter_plot3(mouseID, success_syn_c_RS, cue_syn_c_RS)

success_syn_c_RS_mean = [];
success_syn_c_RS_sem = [];
cue_syn_c_RS_mean = [];
cue_syn_c_RS_sem = [];
for session_num = 1:length(mouseID)
    
%     success_syn_c_RS_mean(session_num) = nanmean(success_syn_c_RS{session_num});
%     success_syn_c_RS_sem(session_num) = nanstd(success_syn_c_RS{session_num})/sqrt(sum(~isnan(success_syn_c_RS{session_num})));
%     cue_syn_c_RS_mean(session_num) = nanmean(cue_syn_c_RS{session_num});
%     cue_syn_c_RS_sem(session_num) = nanstd(cue_syn_c_RS{session_num})/sqrt(sum(~isnan(cue_syn_c_RS{session_num})));

    %removes all std values = 0 
    success_syn_c_RS_mean(session_num) = nanmean(success_syn_c_RS{session_num}(success_syn_c_RS{session_num}>0));
    success_syn_c_RS_sem(session_num) = nanstd(success_syn_c_RS{session_num}(success_syn_c_RS{session_num}>0))/sqrt(sum(~isnan(success_syn_c_RS{session_num}(success_syn_c_RS{session_num}>0))));
    cue_syn_c_RS_mean(session_num) = nanmean(cue_syn_c_RS{session_num}(cue_syn_c_RS{session_num}>0));
    cue_syn_c_RS_sem(session_num) = nanstd(cue_syn_c_RS{session_num}(cue_syn_c_RS{session_num}>0))/sqrt(sum(~isnan(cue_syn_c_RS{session_num}(cue_syn_c_RS{session_num}>0))));
    
    errorbarxy(success_syn_c_RS_mean(session_num), cue_syn_c_RS_mean(session_num), success_syn_c_RS_sem(session_num), cue_syn_c_RS_sem(session_num), {'o', 'k', 'k'});
    hold on;
end

%errorbarxy(success_syn_c_RS_mean, cue_syn_c_RS_mean, success_syn_c_RS_sem, cue_syn_c_RS_sem); hold on;
 scatter(success_syn_c_RS_mean, cue_syn_c_RS_mean, [], 'k');
%xx = [1:10:200];
%plot(xx,xx);