function [stims,blanks] = seq2epochs(Seqposition,nOFF,nON)

nFramePerStim = nOFF + nON;

epochs = {};

for iStim = 1:length(Seqposition);
%     length(Seqposition(iStim).ind)
    for iTrial = 1:length(Seqposition(iStim).ind);
        ind = Seqposition(iStim).ind(iTrial);
        
        blanks{iTrial,iStim}= (ind-1)*nFramePerStim + [1:nOFF];
        stims{iTrial,iStim}= (ind-1)*nFramePerStim + nOFF + [1:nON];
        
    end;
end;

return;




