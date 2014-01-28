goodfits_avg_dF = zeros(5,5,3,4);
dFs= [3.2 1.6 .8 .4];
for idF = 1:length(dFs)
for iArea = 1:3
    nexp = all_fits(iArea).nexp;
    goodfits = [];
    for iexp = 1:nexp
        n = all_fits(iArea).expt(iexp).n(1);
        for iCell = 1:n
        	if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1
                if max(max(all_fits(iArea).expt(iexp).bouton(iCell).plotfit,[],1),[],2)<dFs(idF)
                    goodfits = [goodfits; reshape(all_fits(iArea).expt(iexp).bouton(iCell).plotfit', 1, 25)];
                end
            end
        end
    end
    goodfits_avg_dF(:,:,iArea, idF) = reshape(mean(goodfits,1), 5, 5)';
end
end

figure;
start = 1;
for idF = 1:4
    for iArea = 1:3
    subplot(4,3,start)
    imagesq(goodfits_avg_dF(:,:,iArea,idF))
        if idF == 1
            title(areas(iArea,:));
        end
    colormap(gray)
    start = start+1;
    end
end
subplot(4,3,7)
ylabel('dF<0.8')

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_allareas_SFxTF_dFcutoffs.pdf']);
        print(gcf, '-dpdf', fn_out);