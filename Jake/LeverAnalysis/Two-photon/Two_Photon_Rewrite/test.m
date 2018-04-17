for i = 1:size(success_resp_RS_2,2)
    sResp = success_resp_RS_2{i};
    sResp(isnan(sResp)) = [];
    sHold = succ_hold_dur{i};
    sHold(isnan(sHold)) = [];
    [BS,idx] = histc(sHold,0:0.25:max(sHold)+0.25);
    sHold_bins_all{i} = accumarray(idx(:),sHold,[],@nanmean);
    sResp_std_all{i} = accumarray(idx(:),sHold,[],@sem);
    sResp_bins_all{i} = accumarray(idx(:),sResp,[],@nanmean);
end