function expt = setupstim(expt, off, on, nstim, nframes)
% expt = setupstim(expt,off, on, nstimpertrial, nframestotal)
% getepochs is removed by Kenichi 3/28/08

expt.off = off;
expt.on = on;
expt.stimdur = expt.off+expt.on;
expt.nstim = nstim;
expt.trialdur = expt.nstim*expt.stimdur;
expt.ntrials = floor(nframes/expt.trialdur);
expt.dur = expt.trialdur*expt.ntrials;

%[expt.stims,expt.blanks] = getepochs (expt.off, expt.on, expt.nstim, expt.ntrials);

return;