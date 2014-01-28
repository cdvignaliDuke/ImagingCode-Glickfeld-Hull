%% identify components
newdir_runA = '110324_LG26\run4';
expt_runA = frGetExpt(newdir);
ics_runA = load(fullfile(expt.dirs.analrootpn,expt.filenames.ica));
cs_A = permute(ics_runA.filters,[2,3,1]);
sm_A = stackFilter(cs_A,1.5);

newdir_runB = '110324_LG26\run5';
expt_runB = frGetExpt(newdir);
ics_runB = load(fullfile(expt.dirs.analrootpn,expt.filenames.ica));
cs_B = permute(ics_runB.filters,[2,3,1]);
sm_B = stackFilter(cs_B,1.5);

xc = zeros(20,20);
for ic_A = 1;
    for ic_B = 1;
        xc(ic_A, ic_B) = xcorr(sm_A(ic_A),sm_B(ic_B));
    end
end
