% finds and maps responsive cells for each condition
% 1. calculate variability by trial over base and resp windows
% 2. ttest for significant responses
% 3. define cells by response to event

nCells = size(success_movie,2);

%% 1. calculate response amplitude and variablity by trial over base and resp windows
%by release
base_release_window = 1:pre_release_frames;
resp_release_window = pre_release_frames+1:pre_release_frames+2;

%success
success_base = squeeze(mean(success_movie(:,:,base_release_window),3));
success_resp = squeeze(mean(success_movie(:,:,resp_release_window),3));
success_resp_avg = mean((success_resp-success_base),1);
success_resp_sem = std((success_resp-success_base),[],1)./sqrt(size(success_resp,1));

%fail
fail_base = squeeze(mean(fail_movie(:,:,base_release_window),3));
fail_resp = squeeze(mean(fail_movie(:,:,resp_release_window),3));
fail_resp_avg = mean((fail_resp-fail_base),1);
fail_resp_sem = std((fail_resp-fail_base),[],1)./sqrt(size(fail_resp,1));

%by press
base_press_window = 1:pre_press_frames-5;
resp_press_window = pre_press_frames-1:pre_press_frames;

press_base = squeeze(mean(press_long_movie(:,:,base_press_window),3));
press_resp = squeeze(mean(press_long_movie(:,:,resp_press_window),3));
press_resp_avg = mean((press_resp-press_base),1);
press_resp_sem = std((press_resp-press_base),[],1)./sqrt(size(press_resp,1));

save([dest_sub '_cell_resp.mat'], 'base_release_window', 'resp_release_window', 'base_press_window', 'resp_press_window', 'press_base', 'press_resp', 'success_base', 'success_resp', 'fail_base', 'fail_resp');

%% 2. ttest for significant responses
[success_h, success_p] = ttest(success_base, success_resp, 'dim', 1, 'tail', 'both');
[fail_h, fail_p] = ttest(fail_base, fail_resp, 'dim', 1, 'tail', 'both');
[press_h, press_p] = ttest(press_base, press_resp, 'dim', 1, 'tail', 'both');

%% 3. define cells by response to event
success_resp_cells = find(success_h);
fail_resp_cells = find(fail_h);
press_resp_cells = find(press_h);
fail_only_cells = fail_resp_cells(find(ismember(fail_resp_cells, success_resp_cells)==0));
success_only_cells = success_resp_cells(find(ismember(success_resp_cells, fail_resp_cells)==0));
fail_and_success_cells = find(and(fail_h, success_h));
press_and_success_cells = find(and(success_h,press_h));
press_notsuccess_cells = press_resp_cells(find(ismember(press_resp_cells, success_resp_cells)==0));
success_notpress_cells = success_resp_cells(find(ismember(success_resp_cells, press_resp_cells)==0));
noresponse_cells = find((success_h+fail_h+press_h)==0);
save([dest_sub '_cell_categories.mat'], 'success_resp_cells', 'fail_resp_cells', 'press_resp_cells','noresponse_cells', 'fail_only_cells', 'success_only_cells', 'fail_and_success_cells', 'press_and_success_cells', 'press_notsuccess_cells', 'success_notpress_cells');

%test success vs fail responses (in cells that respond to both)
[fail_v_success_h, fail_v_success_p] = ttest(fail_resp_avg(1,fail_and_success_cells),success_resp_avg(1,fail_and_success_cells), 'tail', 'left');
save([dest_sub '_pvals.mat'], 'fail_v_success_p', 'success_p', 'fail_p', 'press_p');

%% 4. plotting
figure;
tt = ((-pre_release_frames:post_release_frames).*double(ifi))./1000;
errorbar(tt,squeeze(mean(mean(success_movie(:,:,:),1),2)),squeeze(std(mean(success_movie(:,:,:),2),[],1))./sqrt(size(success_movie,1)), '-k');
hold on
errorbar(tt,squeeze(mean(mean(fail_movie(:,:,:),1),2)),squeeze(std(mean(fail_movie(:,:,:),2),[],1))./sqrt(size(fail_movie,1)), '-r');
hold on
errorbar(tt,squeeze(mean(mean(press_movie(:,:,6:end),1),2)),squeeze(std(mean(press_movie(:,:,6:end),2),[],1))./sqrt(size(press_movie,1)), '-c');
title([date ' ' mouse ' Avg resp: Success (black), Fail (red), Press (cyan)']);
print([dest_sub 'Summary_allcell_responses.eps'], '-depsc');
print([dest_sub 'Summary_allcell_responses.pdf'], '-dpdf');

figure;
all_resp = [success_resp_avg fail_resp_avg press_resp_avg];
cmax = max(all_resp,[],2);
cmin = min(all_resp,[],2);
success_mask = mask_final;
sz = size(sm);
for ic = 1:nCells
    if success_h(ic)
        success_mask(find(success_mask==ic))=success_resp_avg(1,ic);
    else
        success_mask(find(success_mask==ic))=0;
    end
end
success_mask = reshape(success_mask,[sz(1) sz(2)]);
subplot(1,3,1)
imagesc(success_mask)
clim([cmin cmax])
title('success')

fail_mask = mask_final;
sz = size(sm);
for ic = 1:nCells
    if fail_h(ic)
        fail_mask(find(fail_mask==ic))=fail_resp_avg(1,ic);
    else
        fail_mask(find(fail_mask==ic))=0;
    end
end
fail_mask = reshape(fail_mask,[sz(1) sz(2)]);
subplot(1,3,2)
imagesc(fail_mask)
clim([cmin cmax])
title('failure')

press_mask = mask_final;
sz = size(sm);
for ic = 1:nCells
    if press_h(ic)
        press_mask(find(press_mask==ic))=press_resp_avg(1,ic);
    else
        press_mask(find(press_mask==ic))=0;
    end
end
press_mask = reshape(press_mask,[sz(1) sz(2)]);
subplot(1,3,3)
imagesc(press_mask)
clim([cmin cmax])
title('press')

suptitle([mouse ' ' date ' cell responses']);
print([dest_sub '_cell_responses_FOV.eps'], '-depsc');
print([dest_sub '_cell_responses_FOV.pdf'], '-dpdf');
