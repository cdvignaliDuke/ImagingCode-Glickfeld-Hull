rshiftdfOvF_run_last1s = zeros(size(pairs,1),size(dfOvF_run_last1s,2)*(size(dfOvF_run_last1s,2)-1));
pshift_dfOvF_run_last1s = zeros(size(pairs,1),size(dfOvF_run_last1s,2)*(size(dfOvF_run_last1s,2)-1));%pairs*trials combinations, combinations = ntrials*ntrials-1 (C ntrial,2)
%shift predictor
for p = 1:size(pairs,1)
    a = 0;
    p
    for t1 = 1:size(dfOvF_run_last1s,2)
        for t2 = 1:size(dfOvF_run_last1s,2)
            if t1~=t2 %for cellA trial1, corr with cellB trial2-N
                a = a+1;
                a
                [rshiftdfOvF_run_last1s(p,a),pshift_dfOvF_run_last1s(p,a)] = corr(dfOvF_run_last1s(:,t1,pairs(p,1)),dfOvF_run_last1s(:,t2,pairs(p,2)));
            end
        end
    end
end
