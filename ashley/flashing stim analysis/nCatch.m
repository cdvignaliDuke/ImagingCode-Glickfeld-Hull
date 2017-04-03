function n = nCatchExpt(expt,mouse)

ncatch = zeros(size(expt));

i = 1;
for im = 1:size(mouse,2)
    for ie = 1:size(mouse(im).expt,2)
        if mouse(im).expt(ie).info.isCatch
            ncatch(i) = length(mouse(im).expt(ie).align(3).av(1).outcome(1).stimResp);
        end
        i = i+1;
    end
end

n = ncatch;
end