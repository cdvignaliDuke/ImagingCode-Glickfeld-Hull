
nMice = length(mice);

expt_tuning = struct;
for im = 1:nMice
    ne = size(mouse(im).expt,2);
    for iexp = 1:ne
        if im == 1 & iexp == 1
            exptN = 1;
        else
            exptN = exptN+1;
        end
        d = mouse(im).expt(iexp);
        expt_tuning(exptN).expt_name = [mouse(im).mouse_name '-' d.date];
        
        
    end
end