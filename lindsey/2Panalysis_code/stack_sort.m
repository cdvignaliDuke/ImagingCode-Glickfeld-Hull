%load sorting data
if length(userun) == 1
    seqfile =[date '_' mouse '_run' num2str(userun) '_Seqposition.mat'];
    load(fullfile(outDir,'analysis',seqfile));
    Big_Seqposition = Seqposition;
else
    seqfile = [date '_' mouse '_run' num2str(userun) '_Big_Seqposition.mat'];
    load(fullfile(outDir,'analysis',seqfile));
end



stack = [];
if P ==1
    for iRun = 1:length(userun);
        substack = readtiff(fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '.tif']));
        stack = cat(3, stack, substack);
    end
elseif P==2
    if nPlanes ==1
        for iRun = 1:length(userun);
            substack = uint16(readtiff(fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '_dec_reg.tif'])));
            stack = cat(3, stack, substack);
        end
    else
        for iRun = 1:length(userun);
            substack = uint16(readtiff(fullfile(outDir, [date '_' mouse '_run' num2str(userun(iRun)) '_plane' num2str(iPlane) '_reg.tif'])));
            stack = cat(3, stack, substack);
        end
    end
end
clear('substack');

% if P ==1 
%     if nOFF>10
%         extra = nOFF-10;
%         nrep = size(stack,3)/(nON+nOFF);
%         ind_mat = zeros(1,nrep*extra);
%         ind = 1;
%         start = 1;
%         for irep = 1:nrep
%             ind_mat(1,ind:ind+extra-1) = start:start+extra-1;
%             start= start+nON+nOFF;
%             ind = ind+extra;
%         end
%         stack(:,:,ind_mat) = [];
%         nOFF = 10;
%     end
% end

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

%resort blanks
nblanks = length(Big_Seqposition(end).ind);
for iblank = 1:nblanks;
    ind = Big_Seqposition(end).ind(iblank);
        if begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)) > size(stack_sorted,3);
            stack_sorted(:,:,start:start-1+((nOFF+nON)/nPlanes)-begin+1) = stack(:,:,begin+((ind-1)*((nOFF+nON)/nPlanes)):((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)));
            stack_sorted(:,:,start+((nOFF+nON)/nPlanes)-begin+1:start+((nOFF+nON)/nPlanes)-1) = stack(:,:,1:begin-1);
        else
        stack_sorted(:,:,start:start-1+((nOFF+nON)/nPlanes)) = stack(:,:,begin+((ind-1)*((nOFF+nON)/nPlanes)):begin-1+((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)));
        end
    start = start+((nOFF+nON)/nPlanes);
end
if start<size(stack_sorted,3);
    stack_sorted(:,:,start:end)=[];
end
%create index of number of reps per stimulus
stim_reps = zeros(1,nCond+1);
for iCond = 1:nCond;
    stim_reps(1,iCond) = length(Big_Seqposition(iCond).ind);
    stim_reps(1,end) = length(Big_Seqposition(end).ind);
end

fn2 = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
save(fn2, 'stim_reps');

%average responses
siz = size(stack_sorted);
stack_dF_allavg = zeros(siz(1), siz(2), nCond+1);
start = 0;

if P ==1;
    if blanks == 1;
        Cond = nCond +1;
    else
        Cond = nCond;
    end
else
    Cond = nCond+1;
end

for iCond = 1:Cond;
    nRep = length(Big_Seqposition(iCond).ind);
    rep_dF = zeros(siz(1), siz(2),nRep);
    for iRep = 1:nRep
        rep_base = mean(stack_sorted(:,:,start+pre_win(1):start+pre_win(2)),3);
        rep_resp = mean(stack_sorted(:,:,start+post_win(1):start+post_win(2)),3);
        rep_dF(:,:,iRep) = (rep_resp-rep_base)./rep_base;
        start = ((nOFF+nON)/nPlanes)+start;
    end
    stack_dF_allavg(:,:,iCond) = mean(rep_dF,3);
end

if nPlanes == 1
fn_out = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
fn = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted.tif']);
else
fn_out = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_plane' num2str(iPlane) '_stack_dF_all.tif']);
fn = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_plane' num2str(iPlane) '_sorted.tif']);
end
writetiff(stack_dF_allavg, fn_out);
writetiff(stack_sorted,fn);

if P==1
    if TFSFetc == [1:2];
        if nCond== 9
            low = [1 2 4];
            med = [3 5 7];
            hi = [6 8 9];
        elseif nCond == 25
            low = [1 2 3 6 7 11];
            med = [4 5 8 9 10 12 13 14 16 17 18 21 22];
            hi = [15 19 20 23 24 25];
        end


        avg_low = mean(stack_dF_allavg(:,:,low),3);
        avg_med = mean(stack_dF_allavg(:,:,med),3);
        avg_hi = mean(stack_dF_allavg(:,:,hi),3);
        avg_lowmedhigh = cat(3, avg_low, avg_med, avg_hi);

        fn_out = fullfile(outDir, 'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_lowmedhigh.tif']);
        writetiff(avg_lowmedhigh, fn_out);

%extract variables
        if P==1;
            if nPlanes == 1 | iPlane == 1;
                clear vars
                clear uvars
                for iS = 1:length(Big_Seqposition)
                    for iVar = TFSFetc;
                        vars(iS,find(iVar==TFSFetc))=Big_Seqposition(iS).TFSFetc(iVar);
                    end
                end

                if size(vars,2)== 1
                    uvars = unique(vars);
                    uvar = uvars;
                else
                    if size(unique(vars(:,1))) == size(unique(vars(:,2)))
                        for iVar = 1:2;
                            uvars(:,iVar) = unique(vars(:,iVar));
                        end
                        nvars = size(uvars);
                        uvars(find(uvars==0)) = [];
                        uvar = reshape(uvars,[nvars(1)-1 nvars(2)]);
                    end
                end
                fn3 = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_vars.mat']);
                save(fn3, 'uvar');
            end
        end
    end
end

clear('stack_avg');
% clear('stack_sorted');
clear('stack');
