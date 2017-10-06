function normTCEaCycEaCell = getNormTCByCondition(tcEaCycEaCell,normCycle,...
    normCondition,respwin)

% [ncycles,nconditions] = size(tcEaCycEaCell);
% 
% normTCEaCycEaCell = cell(ncycles,nconditions);
% for i = 1:nconditions
%     normTCEaCycEaCell(:,i) = cellfun(@(x) ...
%         x./mean(tcEaCycEaCell{normCycle,i}(respwin,:),1),...
%         tcEaCycEaCell(:,i),'unif',0);
% end



normTCEaCycEaCell = cell(size(tcEaCycEaCell));
for iav = 1:2
    normTCEaCycEaCell(:,iav) = cellfun(@(x)...
        x./mean(tcEaCycEaCell{normCycle,iav}(respwin,:),1),...
        tcEaCycEaCell(:,iav),'unif',0);
end

% normTCEaCycEaCell = cellfun(@(x)...
%     x./mean(tcEaCycEaCell{normCycle,normCondition}(respwin,:),1),...
%     tcEaCycEaCell,'unif',0);

end