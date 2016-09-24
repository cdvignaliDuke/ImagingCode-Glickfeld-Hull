Nshuf = 500;
SFs = cell2mat(input.tGratingSpatialFreqCPD);
TFs = cell2mat(input.tGratingTemporalFreqCPS);
SF_vec0 = flipud(unique(SFs)); %flipped to have low to high SF in square  %flipud
TF_vec0 = unique(TFs);
nSF = length(SF_vec0);
nTF = length(TF_vec0);

[tftf,sfsf] = meshgrid(TF_vec0,SF_vec0); 
grid2.sfsf = sfsf;
grid2.tftf = tftf;

dSF = median(diff(log2(SF_vec0)));
dTF = median(diff(log2(TF_vec0)));
SF_vec00 = log2(SF_vec0(1)):(dSF/10):log2(SF_vec0(end));
TF_vec00 = log2(TF_vec0(1)):(dTF/10):log2(TF_vec0(end));
[sfsf00,tftf00]=meshgrid(SF_vec00,TF_vec00);
grid2.sfsf00 = sfsf00;
grid2.tftf00 = tftf00;

Ind_struct = [];
iCond = 1;
for iSF = 1:nSF
    ind1 = find(SFs == SF_vec0(iSF));
    for iTF = 1:nTF
        ind2 = find(TFs == TF_vec0(iTF));
        Ind_struct(iCond).all_trials = intersect(ind1,ind2);
        iCond = iCond+1;
    end
end
    
    
Fit_struct = [];
    for count_shuf = 0:Nshuf
        fprintf('.')
        Im_mat_USE = zeros(size(resp_dF,1), 25);
        Im_mat_std = zeros(size(resp_dF,1), 25);
        for iCond = 1:25        
            ind_all = Ind_struct(iCond).all_trials;
            if count_shuf > 0 %resample with replacement, don't resample by trial for now because running-rejection may be uneven for various trials..
                ind_all_1 = ind_all(randsample(length(ind_all),length(ind_all),1));
            else
                ind_all_1 = ind_all;        
            end
            Im_mat_USE(:,iCond) = mean(resp_dF(:,ind_all_1),2);
        end

        start = 1;
        for iCell = 1:size(resp_dF,1);
            a = Im_mat_USE(iCell,:);
            b = reshape(a',length(SF_vec0),length(TF_vec0));
            %b2 = b( ind_SFuse(:,1),ind_TFuse(:,1));
            data = b';
            ind0 = find(data<0);
            data(ind0) = NaN;
            if count_shuf == 0
                PLOTIT_FIT = 0;
                SAVEALLDATA = 1;
                Fit_2Dellipse_LG
                eval(['Fit_struct(iCell).True.s_',' = s;']);
            else
                SAVEALLDATA = 0;
                PLOTIT_FIT = 0;
                Fit_2Dellipse_LG
                eval(['Fit_struct(iCell).Shuf(count_shuf).s_',' = s;']);
            end
        end
    end

    fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct.mat']);   
    save(fn_out, 'Fit_struct')


if Nshuf>1;
    for iCell = 1:n_pix
        if ~isempty(Fit_struct(iCell).True)                
            eval(['tmp = Fit_struct(iCell).True.s_.x;']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.SFhicut_50];']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.TFhicut_50];']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.SFhicut_10];']);
            eval(['tmp = [tmp Fit_struct(iCell).True.s_.TFhicut_10];']);
            fit_true_vec(iCell,:) = tmp;
        end
    end

    for count_shuf = 1:Nshuf
        for iCell = 1:n_pix
            if ~isempty(Fit_struct(iCell).Shuf)
                eval(['tmp = Fit_struct(iCell).Shuf(count_shuf).s_.x;']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.SFhicut_50];']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.TFhicut_50];']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.SFhicut_10];']);
                eval(['tmp = [tmp Fit_struct(iCell).Shuf(count_shuf).s_.TFhicut_10];']);
                %fit is: %A sigma_SF sigma_TF sf0 tf0 xi
                fit_shuf_vec(iCell,:,count_shuf) = tmp;
            end
        end
    end

    Npars = size(fit_shuf_vec,2);
    lbub_fits = zeros(n_pix,Npars,5);
    alpha_bound = .025;
    for iCell = 1:n_pix
        for count2 = 1:Npars
            tmp = squeeze(fit_shuf_vec(iCell,count2,:));
            [i,j] = sort(tmp);
            ind_shuf_lb = ceil(Nshuf*alpha_bound);
            ind_shuf_ub = ceil(Nshuf*(1-alpha_bound));
            lbub_fits(iCell,count2,1) = i(ind_shuf_lb);
            lbub_fits(iCell,count2,2) = i(ind_shuf_ub);
            lbub_fits(iCell,count2,3) = mean(i); 
            lbub_fits(iCell,count2,5) = std(i);
        end
        %now take means from truedata fit:
        lbub_fits(iCell,:,4) = fit_true_vec(iCell,:);
    end
end

lbub_diff = lbub_fits(:,:,2)-lbub_fits(:,:,1);

goodfit = [];
for iCell = 1:n_pix
    if lbub_diff(iCell,4)<2 
        if lbub_diff(iCell,5)<2
            goodfit = [goodfit iCell];
        end
    end
end

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);   
save(fn_out, 'lbub_fits', 'lbub_diff', 'goodfit')