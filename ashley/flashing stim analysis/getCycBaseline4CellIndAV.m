function [avgBaselineEaCycAV,semBaselineEaCycAV] = getCycBaseline4CellIndAV(...
    tcAVEaCellCmlvCyc,cellInd,normCondition,basewin,respwin,cycLengthFr,doResp)
% doResp: 0 = keep all significant cells, 1 = keep only cells with big
% enough response, 2 = keep only cells with big enough suppression

    [ncycles,nexp,nav] = size(tcAVEaCellCmlvCyc);
    avgBaselineEaCycAV = zeros(ncycles,nav);
    semBaselineEaCycAV = zeros(ncycles,nav);
    tcEaCycAV = cell(ncycles,nav);
    for icyc = 1:ncycles
        extraFrames = cycLengthFr*(icyc-1);
        for iexp = 1:nexp
            for iav = 1:nav
                ind = logical(cellInd{iexp,1});
                tcEaCycAV{icyc,iav} = [tcEaCycAV{icyc,iav} ...
                    tcAVEaCellCmlvCyc{icyc,iexp,iav}(:,ind)];
            end 
        end 
        if icyc == 1
            response2firstStim = mean(...
                tcEaCycAV{1,normCondition}(respwin,:),1);
            if doResp == 0
                bigind = true(1,length(response2firstStim));
            elseif doResp == 1
                bigind = response2firstStim > 0.002;
            elseif doResp == 2
                bigind = response2firstStim < 0.002;
            end
        end
        for iav = 1:nav
            avgBaselineEaCycAV(icyc,iav) = mean(mean(...
                tcEaCycAV{icyc,iav}(basewin+extraFrames,bigind),1),2);
            semBaselineEaCycAV(icyc,iav) = ste(mean(...
                tcEaCycAV{icyc,iav}(basewin+extraFrames,bigind),1),2);
        end
    end
end