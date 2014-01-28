areas = ['PM'; 'LM'; 'AL'];
for iArea = 1:3;
P = 2;
matrix = 'SF5xTF5';
image = areas(iArea,:);
inj = 'V1';
nON=12;
nOFF=12;
Nshuf = 500;
SF_vec0 = [.32 .16 .08 .04 .02]; %flipped to have low to high SF in square  %flipud
TF_vec0 = [1 2 4 8 15];

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

sum_base = 'G:\users\lindsey\analysisLG\experiments';
base = 'G:\users\lindsey\analysisLG\active mice';

list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
load(list_fn);
nexp = size(exp_list.mouse_mat,2);
mouse_list = [];
mice = exp_list.mouse_mat;

for iexp = 1:nexp
    mouse = char(exp_list.mouse_mat{iexp});
    date = char(exp_list.date_mat{iexp});
    userun = exp_list.runs_mat{iexp};
    count_prot = exp_list.prot_mat{iexp};
    run = exp_list.run_mat{iexp};
    blanks = exp_list.blanks_mat{iexp};
    dir = exp_list.dir_mat{iexp};
    
    base = 'G:\users\lindsey\analysisLG\active mice';    
    outDir = fullfile(base, mouse,date);

        fn_resp = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp.mat']);
        load(fn_resp);

        fn_mask = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
        load(fn_mask);

        fn_reps = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
        load(fn_reps);

    if dir == 2;
        start = 1;
        for iCond = 1:2:50;
            stim_reps_temp(1,start) = sum(stim_reps(1,iCond:1+iCond),2);
            start = start+1;
        end
        stim_reps = stim_reps_temp;
    end
        
    Ind_struct = [];
    start = 1;
    for iCond = 1:25
        nRep = stim_reps(iCond);
        Ind_struct(iCond).all_trials = [start:start-1+nRep];
        start = start+nRep;
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
    end
end

figure;
start = 1;
fig = 1;
for iCell = 1:n_pix    
    if start>64
        fn_out = fullfile('\\Zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\120112', [date '_' mouse '_run' num2str(userun) '_Data_vs_Fit_' num2str(fig) '.pdf']);
        print(gcf, '-dpdf', fn_out);
        figure;
        fig = 1+fig;
        start = 1;
    end
    subplot(8, 8, start);
    imagesq(Fit_struct(iCell).True.s_.data);
    subplot(8, 8, start+1);
    imagesq(Fit_struct(iCell).True.s_.k2b_plot);
    if lbub_diff(iCell,4)<2 
        if lbub_diff(iCell,5)<2
            title('**');
        end
    end
    colormap(gray);
    start= start+2;
end
fn_out = fullfile('\\Zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\120112', [date '_' mouse '_run' num2str(userun) '_Data_vs_Fit_' num2str(fig) '.pdf']);
print(gcf, '-dpdf', fn_out);

figure;
for iCell = 1:n_pix
    if lbub_diff(iCell,4)<2 
        if lbub_diff(iCell,5)<2
            loglog(2.^lbub_fits(iCell,5,4), 2.^lbub_fits(iCell,4,4),'.k')
            hold on
        end
    end
end
xlim([1 15])
ylim([0.02 0.32])
xlabel('Temporal Frequency')
ylabel('Spatial Frequency')
fn_out = fullfile('\\Zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\120111', [date '_' mouse '_run' num2str(userun) '_SF_vs_TF_scatter.pdf']);
        print(gcf, '-dpdf', fn_out);



pass = find(lbub_diff(:,4)<2 & lbub_diff(:,5)<2);

speed = (2.^lbub_fits(:,5,4))./(2.^lbub_fits(:,4,4));
dF = squeeze(lbub_fits(:,1,4));
dF_pass = dF(pass);
speed_pass = speed(pass);
speed_vec0 = TF_vec0./SF_vec0;

speed_avg = zeros(5,2);
dF_avg = zeros(5,2);
n = zeros(5,1);
for ispeed = 1:5
    if ispeed == 1
        ind = find(speed_pass<=speed_vec0(1,ispeed)+0.00001);
    else
        ind = find(speed_pass<=speed_vec0(1,ispeed)+0.00001 & speed_pass>speed_vec0(1,ispeed-1)+0.00001);
    end
    speed_avg(ispeed,1) = mean(speed_pass(ind,:),1);
    speed_avg(ispeed,2) = std(speed_pass(ind,:),1)./sqrt(size(ind,1));
    dF_avg(ispeed,1) = mean(dF_pass(ind,:),1);
    dF_avg(ispeed,2) = std(dF_pass(ind,:),1)./sqrt(size(ind,1));
    n(ispeed,1) = length(ind);
end
figure;
ploterr(speed_avg(:,1), dF_avg(:,1), speed_avg(:,2), dF_avg(:,2), 'logx')
hold on
for ispeed = 1:5
    text(speed_avg(ispeed,1)-.1*speed_avg(ispeed,1), 0.1, num2str(n(ispeed, 1)));
end
ylim([0 2]);
xlim([0 750]);
xlabel('Speed');
ylabel('dF/F');

fn_out = fullfile('\\Zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\120111', [date '_' mouse '_run' num2str(userun) '_Speed_vs_dF_errorbar.pdf']);
        print(gcf, '-dpdf', fn_out);

figure;
for iCell = 1:n_pix
    if lbub_diff(iCell,4)<2 
        if lbub_diff(iCell,5)<2
            semilogx(speed(iCell,:), lbub_fits(iCell,1,4),'.k')
            hold on
        end
    end
end
xlim([0 1000])
xlabel('Speed')
ylabel('dF/F')

fn_out = fullfile('\\Zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\120111', [date '_' mouse '_run' num2str(userun) '_Speed_vs_dF_scatter.pdf']);
        print(gcf, '-dpdf', fn_out);
