%load sorting data
seqfile = [date '_' mouse '_run' num2str(userun) '_Big_Seqposition.mat'];
load(fullfile(outDir,'analysis',seqfile));

stack = [];
for iRun = 1:length(userun);
    substack = readtiff(fullfile(base,mouse,date, [date '_' mouse '_run' num2str(userun(iRun)) '.tif']));
    stack = cat(3, stack, substack);
end
clear('substack');

stack_sorted = zeros(size(stack) ,'uint16');

%resort stimuli
start = 1;
for iCond = 1:nCond;
    nRep = length(Big_Seqposition(iCond).ind);
    for iRep = 1:nRep;
        ind = Big_Seqposition(iCond).ind(iRep);
        if begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)) > size(stack_sorted,3);
            stack_sorted(:,:,start:start-1+((nOFF+nON)/nPlanes)-begin+1) = stack(:,:,begin+((ind-1)*((nOFF+nON)/nPlanes)):((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)));
            stack_sorted(:,:,start+((nOFF+nON)/nPlanes)-begin+1:start+((nOFF+nON)/nPlanes)-1) = stack(:,:,1:begin-1);
        else
        stack_sorted(:,:,start:start-1+((nOFF+nON)/nPlanes)) = stack(:,:,begin+((ind-1)*((nOFF+nON)/nPlanes)):begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)));
        end
        start = start+((nOFF+nON)/nPlanes);
    end
end

% %resort blanks
% nblanks = length(Big_Seqposition(end).ind);
% for iblank = 1:nblanks;
%     ind = Big_Seqposition(end).ind(iblank);
%         if begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)) > size(stack_sorted,3);
%             stack_sorted(:,:,start:start-1+((nOFF+nON)/nPlanes)-begin+1) = stack(:,:,begin+((ind-1)*((nOFF+nON)/nPlanes)):((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)));
%             stack_sorted(:,:,start+((nOFF+nON)/nPlanes)-begin+1:start+((nOFF+nON)/nPlanes)-1) = stack(:,:,1:begin-1);
%         else
%         stack_sorted(:,:,start:start-1+((nOFF+nON)/nPlanes)) = stack(:,:,begin+((ind-1)*((nOFF+nON)/nPlanes)):begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)));
%         end
%     start = start+((nOFF+nON)/nPlanes);
% end

%create index of number of reps per stimulus
stim_reps = zeros(1,nCond);
for iCond = 1:nCond;
    stim_reps(1,iCond) = length(Big_Seqposition(iCond).ind);
end

%average responses
siz = size(stack_sorted);
stack_avg = zeros(siz(1), siz(2), nON+nOFF, nCond);
stack_avg_base = zeros(siz(1), siz(2), nCond);
stack_avg_resp = zeros(siz(1), siz(2), nCond);
start = 0;
for iCond = 1:nCond;
    for iFrame = 1:(nOFF+nON)
    stack_avg(:,:,iFrame,iCond) = mean(stack_sorted(:,:,iFrame+start:nOFF+nON:start+((nOFF+nON)*stim_reps(iCond))),3);
    end
    stack_avg_base(:,:,iCond) = mean(stack_avg(:,:,pre_win(1):pre_win(2),iCond),3);
    stack_avg_resp(:,:,iCond) = mean(stack_avg(:,:,post_win(1):post_win(2),iCond),3);
    stack_dF_all(:,:,iCond) = (stack_avg_resp(:,:,iCond)-stack_avg_base(:,:,iCond))./stack_avg_base(:,:,iCond);
    start = (nOFF+nON)*stim_reps(iCond)+start;
end

stack_dF = bsxfun(@minus, reshape(stack_avg,[siz(1) siz(2) (nON+nOFF)*(nCond)]), mean(stack_sorted,3));
fn_out = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
writetiff(stack_dF_all, fn_out);

%extract variables
clear vars
clear uvars
for iS = 1:length(Big_Seqposition)
    for iVar = TFSFetc;
        vars(iS,find(iVar==TFSFetc))=Big_Seqposition(iS).TFSFetc(iVar);
    end
end

for iVar = 1:length(TFSFetc);
    uvars(:,iVar) = unique(vars(:,iVar));
end

if TFSFetc <= 3
    nvars = size(uvars);
    uvars(find(uvars==0)) = [];
    uvar = reshape(uvars,[nvars(1) nvars(2)]);
end
if TFSFetc > 3
    uvar = uvars;
end

%save 
fn = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted.tif']);
writetiff(stack_sorted,fn);

fn2 = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
save(fn2, 'stim_reps');

fn3 = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_vars.mat']);
save(fn3, 'uvar');

clear('stack_avg');
clear('stack_sorted');
clear('stack');
