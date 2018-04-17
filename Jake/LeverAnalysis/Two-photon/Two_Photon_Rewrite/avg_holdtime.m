function [sHold_mean, sResp_mean,sResp_std, idx1, idx2] = avg_holdtime(success_resp_RS_2, succ_hold_dur)
for i = 1:size(success_resp_RS_2,2)
    sResp = success_resp_RS_2{i};
    sResp(isnan(sResp)) = [];
    sHold = succ_hold_dur{i};
    sHold(isnan(sHold)) = [];
    [BS,idx] = histc(sHold,0:0.25:6);
    sHold_bins_all{i} = accumarray(idx(:),sHold,[],@nanmean);
    sResp_std_all{i} = accumarray(idx(:),sHold,[],@sem);
    sResp_bins_all{i} = accumarray(idx(:),sResp,[],@nanmean);
    if ~isempty(idx)
        idx_min{i} = min(idx);
        idx_max{i} = max(idx);
    else
        idx_min{i} = [];
        idx_max{i} = [];
    end
end
lv = max(cellfun(@length,sResp_bins_all));
za = zeros(lv,1);
for idx = 1:length(sResp_bins_all)
    cellVal = sResp_bins_all{idx};
  if length(cellVal) < lv
      sResp_bins_all{idx} = za;
      sResp_bins_all{idx}(1:length(cellVal)) = cellVal;
  end
end

lv = max(cellfun(@length,sHold_bins_all));
za = zeros(lv,1);
for idx = 1:length(sHold_bins_all)
    cellVal = sHold_bins_all{idx};
  if length(cellVal) < lv
      sHold_bins_all{idx} = za;
      sHold_bins_all{idx}(1:length(cellVal)) = cellVal;
  end
end

sResp = cell2mat(sResp_bins_all);
sHold = cell2mat(sHold_bins_all);
sHold_mean = sum(sHold, 2)./sum(sHold~=0, 2);
sHold_mean(isnan(sHold_mean)) = 0;
sResp_mean = sum(sResp, 2)./sum(sResp~=0, 2);
sResp_mean(isnan(sResp_mean)) = 0;
sResp_std = std(sResp, [], 2)./sqrt(size(sResp,2));
sResp_std(isnan(sResp_std)) = 0;
idx1 = min(cell2mat(cellfun(@int64,idx_min,'UniformOutput',0)));
idx2 = max(cell2mat(idx_max));
end