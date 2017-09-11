function [tcEaCycAV,avgResponseEaCycAV,semResponseEaCycAV] = ...
    getCycResponse4CellIndAV(tcAVEaCell,cellInd,normCondition,respwin,doResp)
% doResp: 0 = keep all significant cells, 1 = keep only cells with big
% enough response, 2 = keep only cells with big enough suppression
    [ncycles,nexp,nav] = size(tcAVEaCell);
    avgResponseEaCycAV = zeros(ncycles,nav);
    semResponseEaCycAV = zeros(ncycles,nav);
    tcEaCycAV = cell(ncycles,nav);
    for icyc = 1:ncycles
        for iexp = 1:nexp
            for iav = 1:nav
                ind = logical(cellInd{iexp,1});
                tcEaCycAV{icyc,iav} = [tcEaCycAV{icyc,iav} ...
                    tcAVEaCell{icyc,iexp,iav}(:,ind)];
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
            avgResponseEaCycAV(icyc,iav) = mean(mean(...
                tcEaCycAV{icyc,iav}(respwin,bigind),1),2);
            semResponseEaCycAV(icyc,iav) = ste(mean(...
                tcEaCycAV{icyc,iav}(respwin,bigind),1),2);
        end
    end
    tcEaCycAV = cellfun(@(x) x(:,bigind),tcEaCycAV,'unif',0);
end