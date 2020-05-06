%% compare size-tuning across visual areas
% script that loads given experiments and sorts by visual area to compare
% size tuning parameters across V1 + HVAs

% given csv .txt file for experiment list
% includes: date, mouse, run_str(size-tuning exp), HVA, indicator

% data to load:
% ?ret results - lbub_fits (w/ RF location), goodfit cells (maybe?)
% ?size tuning run input - stim location (maybe?)
% _sizeTuneData.mat - sizeTune, sizeMean, sizeSEM, cellDists
% _sizeFitResults_SP.mat - sizeFits struct (all cons)
% _Fit_struct.mat - Fit_struct (highest con, with shuffles)
% _lbub_fits - lbub_fits, goodfit_ind_size

%% load csv to define exp

clear all; clc;
fn_in = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';
fn_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Syt7';
expfile = fullfile(fn_out, 'syt7Ret_experiment_list.txt');
fID = fopen(expfile);
head = textscan(fID,'%s%s%s%s%s%s',1,'delimiter',',');
head = vertcat(head{:});
temp = textscan(fID,'%s%s%s%s%s%s','delimiter',',','HeaderLines',1);
temp = horzcat(temp{:});
expdata = cell2table(temp,'VariableNames',head);
nExp = size(expdata,1);
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);
fprintf(['Retinotopy Syt7 comparison analysis - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

%% load each experiment and concatenate data

fprintf('\nBegin loading and concatentating experiment data...\n')
% sizeTuneData
fprintf('Loading retinotopy data\n')
lbub_fits_all = [];
goodfit_ind_all = [];

t = 0;
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    ret_run_str = expdata.ret_run_str{i};
    filename = fullfile(fn_in, [date '_i' mouse], [ret_run_str], [mouse '_' date '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' ret_run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind')
    
    
    nCellsExp(i) = size(lbub_fits,1);
    lbub_fits_all = cat(1,lbub_fits_all,lbub_fits);
    tempinds = sum(nCellsExp(1:i-1)) + goodfit_ind; % offset by # cells in previous exps
    goodfit_ind_all = [goodfit_ind_all tempinds];
    
    fprintf('done\n')
end

nCellsTot = sum(nCellsExp);
fprintf('%d cells loaded (%d goodfit)\n',nCellsTot,length(goodfit_ind_all))
expInd = [];
expArea = [];
expGenotype = [];
for i = 1:nExp
    expInd = [expInd repmat(i,1,nCellsExp(i))];
    expArea = [expArea repmat(expdata(i,4).area, 1,nCellsExp(i))];
    expGenotype = [expGenotype repmat(expdata(i,6).genotype, 1,nCellsExp(i))];
end


%% compare genotypes
% making figs 1-7
genotypes = unique(expGenotype);
nGenotype = length(genotypes);

nExp_genotype = zeros(1,nGenotype);
nCells_genotype = zeros(1,nGenotype);

close all
% choose figs: 1=modelcounts; 2=averagecurves; 3=prefSize; 4=suppInd; 5=conresp; 6=ex.cells; 7=medianfits;
legStrs = strings(1,nGenotype); legStrs2=legStrs;
for i = 1:nGenotype
    fprintf(['Genotype #' num2str(i) ' : ' char(genotypes(i)) '\n'])
    % select exps matching area
    expIndi = find(cellfun(@(x) strcmp(x,genotypes(i)), expdata.genotype, 'UniformOutput', 1));
    % find cells with correct exp inds, take only good fit cells
    ind = intersect(find(ismember(expInd,expIndi)),goodfit_ind_all);
    
    
    nExpi = length(expIndi);
    nCellsi = length(ind);
    nExp_genotype(i) = nExpi;
    lbub_fits = lbub_fits_all(ind,:,:); %cell,par,val (low up mean true stdev)
    sigmax = lbub_fits(:,2,4);
    sigmay = lbub_fits(:,3,4);
    RFsize_all = 2*sqrt(2*log(2))*geomean([sigmax sigmay],2);
    
    legStrs(i)=sprintf('%s (n=%d, n_{exp}=%d)',cell2mat(genotypes(i)),nCellsi,nExpi);
    
    figure(1);if i==1;clf;end %figure 1 = proportions of model2 by con
    subplot(1, nGenotype, i)
    scatter(ones(size(RFsize_all)), RFsize_all, 'o')
    hold on
    errorbar(1,mean(RFsize_all,1), std(RFsize_all,[],1)./sqrt(nCellsi),'o');
    title({sprintf('Genotype:%s',cell2mat(genotypes(i)));['(n=' num2str(nCellsi) ', n_{exp}=' num2str(nExpi) ')']})
    ylabel('FWHM (deg)')
    ylim([0 30])
end
figure(1)
print(fullfile(fn_out, 'syt7_RFsize.pdf'),'-dpdf','-bestfit')