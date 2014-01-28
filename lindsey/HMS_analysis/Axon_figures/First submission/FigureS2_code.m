fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 'S2';
areas = strvcat('PM', 'LM', 'AL');
area_order = [2;3;1];
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);


%example mouse
mouse ='M14';
date = '111127';
userun = 1:4;
iexp = 12;
iArea = 1;
dirs = 2;

base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse,date);

%example dF/F for speeds

fn_df =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_stack_dF_np.tif']);
stack = readtiff(fn_df);

siz = size(stack);
dF_nodir = zeros(siz(1), siz(2), 25);
start = 1;
for iCond = 1:25
    dF_nodir(:,:,iCond) = mean(stack(:,:,start:start+1),3);
    start = start+2;
end

figure;
dF_max = max(max(max(dF_nodir,[],1),[],2),[],3);
for iframe = 1:25
    subplot(5,5,iframe)
    imagesq(dF_nodir(:,:,iframe)./dF_max);
    colormap(gray)
    caxis([0 .05])
end

fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_all_areas_' mouse '_' date '_all_dF_resp.ps']);
        print(gcf, '-depsc2', fn_out);

%example correlation matrix
x = 1:5;
y = 5:-1:1;
[x_grid y_grid] = meshgrid(x,y);
x_grid_long = reshape(x_grid', [25 1]);
y_grid_long = reshape(y_grid', [25 1]);

fn_resp = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp.mat']);
load(fn_resp);
fn_reps = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
load(fn_reps);
fn_lbub = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
load(fn_lbub);


n2 = all_fits(1).expt(iexp).n(2);
r_real = zeros(n2);
for ipix1 = 1:n2
    for ipix2 = 1:n2
        r_real(ipix1,ipix2)= triu2vec(corrcoef(resp_dFoverF(goodfit_ind(:,ipix1),:),resp_dFoverF(goodfit_ind(:,ipix2),:)));
    end
end
figure; 
subplot(2,2,1); 
imagesq(r_real);
[order, groups] = corrsort(r_real);
r_order = r_real(order,order);
colorbar
subplot(2,2,2);
imagesq(r_order);
colorbar

fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_all_areas_' mouse '_' date '_corr_matrix.ps']);
        print(gcf, '-depsc2', fn_out);

n2 = all_fits(1).expt(iexp).n(2);
r_strip = zeros(n-1,1);
for ipix = 1:n-1
    r_strip(ipix,:) = r_order(ipix,ipix+1);
end

pairs = find(r_strip>.75);

n = all_fits(1).expt(iexp).n(1);
plotfit = [];
speed = [];
for iCell = 1:n
    if all_fits(1).expt(iexp).bouton(iCell).goodfit == 1
        plotfit_reshape =  reshape(all_fits(1).expt(iexp).bouton(iCell).plotfit', [25 1]);
        plotfit_norm = (plotfit_reshape./max(plotfit_reshape,[],1))*20;
        plotfit = [plotfit plotfit_norm];
        speed = [speed all_fits(1).expt(iexp).bouton(iCell).speed];
    end
end


figure;
start =1;
fignum = 1;
for ipair = 1:length(pairs)
    cellA = order(:,pairs(ipair));
    subplot(6,6,start)
    for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), plotfit(iCond,cellA).^2,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
    end
    title(num2str(pairs(ipair)));
    start = start+1;
    if start==37
        fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_all_areas_' mouse '_' date '_boutongroups_Figure' num2str(fignum) '.ps']);
        print(gcf, '-depsc2', fn_out);
        figure;
        start = 1;
        fignum = fignum+1;
    end
    if pairs(ipair+1)>(pairs(ipair)+1)
         cellB = order(:,pairs(ipair)+1);
         subplot(6,6,start)
         for iCond = 1:25
            scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), plotfit(iCond,cellB).^2,'.k')
            hold on
            axis square
            axis off
            xlim([0 6])
            ylim([0 6])
         end
         title(num2str(pairs(ipair)+1));
         start = start+1;
    elseif pairs(ipair) == pairs(end)
         cellB = order(:,pairs(ipair)+1);
         subplot(6,6,start)
         for iCond = 1:25
            scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), plotfit(iCond,cellB).^2,'.k')
            hold on
            axis square
            axis off
            xlim([0 6])
            ylim([0 6])
         end
         title(num2str(pairs(ipair)+1));
         start = start+1;
    end
    if start==37
        fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_all_areas_' mouse '_' date '_boutongroups_Figure' num2str(fignum) '.ps']);
        print(gcf, '-depsc2', fn_out);
        figure;
        start = 1;
        fignum = fignum+1;
    end
end
x = 37-start;
for ix = 1:x
    subplot(6,6,start)
     for iCond = 1:25
        scatter(x_grid_long(iCond,1),y_grid_long(iCond,1), 1000,'.k')
        hold on
        axis square
        axis off
        xlim([0 6])
        ylim([0 6])
     end
     start = start+1;
end
fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_all_areas_' mouse '_' date '_boutongroups_Figure' num2str(fignum) '.ps']);
        print(gcf, '-depsc2', fn_out);
        
%single trial time courses for individual boutons
fn_reg = fullfile(outDir,[date '_' mouse '_run1_dec_reg.tif']);
reg = readtiff(fn_reg);

TC = [];
pos = [];
start= 1;
for ipair = 1:length(pairs)
    cell = order(:,pairs(ipair));
    iCell = goodfit_ind(:,cell);
    pos(:,start) = all_fits(iArea).expt(iexp).bouton(iCell).pos;
    suby = pos(1,start)-1:pos(1,start)+1;
    subx = pos(2,start)-1:pos(2,start)+1;
    TC(:,start) = squeeze(mean(mean(reg(suby,subx,:),2),1));
    start = start+1;
    if ipair == length(pairs)
        cell = order(:,pairs(ipair)+1);
        iCell = goodfit_ind(:,cell);
        pos(:,start) = all_fits(iArea).expt(iexp).bouton(iCell).pos;
        suby = pos(1,start)-1:pos(1,start)+1;
        subx = pos(2,start)-1:pos(2,start)+1;
        TC(:,start) = squeeze(mean(mean(reg(suby,subx,:),2),1));
        start = start+1;
    elseif pairs(ipair+1)>(pairs(ipair)+1)
        cell = order(:,pairs(ipair)+1);
        iCell = goodfit_ind(:,cell);
        pos(:,start) = all_fits(iArea).expt(iexp).bouton(iCell).pos;
        suby = pos(1,start)-1:pos(1,start)+1;
        subx = pos(2,start)-1:pos(2,start)+1;
        TC(:,start) = squeeze(mean(mean(reg(suby,subx,:),2),1));
        start = start+1;
    end   
end

TC_dF = bsxfun(@minus, TC,mean(TC,1));
TC_dFoverF = bsxfun(@rdivide, TC_dF, mean(TC,1));

seqfile =[date '_' mouse '_run1_Seqposition.mat'];
    load(fullfile(outDir,'analysis',seqfile));
stim_time = [13:24:size(reg,3)];
stim_times = stim_time;
stim_times(:,Seqposition(51).ind)=[];
figure; tcOffsetplot(TC_dFoverF(:,[1 2 4 5 13 14 15 16]));
figure; tcOffsetplot(TC_dFoverF(:,[25 26 27 28 31 32 33 34]));
figure; tcOffsetplot(TC_dFoverF(:,[37 38 39 40 46 47 53 54]));
figure; tcOffsetplot(TC_dFoverF(:,[56 57 69 70 71 72 73]));
hold on
for itime = 1:length(stim_times)
    x = stim_time(:,itime);
    y = -2:1:20;
    plot(x*ones(length(y)),y,':k');
    hold on
end
xlim([1370 2339])
ylim([-2 18]) 
fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_all_areas_' mouse '_' date '_boutonTCs_row3.ps']);
        print(gcf, '-depsc2', fn_out);

pair_list = [1 3 11 12 21 22 25 26 29 30 35 40 43 53 54 55];
figure; 
plot(r_strip);
xlim([1 length(r_strip)]);
hold on
plot(1:1:length(r_strip), .75, ':k')
for ipair = 1:16
    x = pairs(pair_list(ipair));
    y = 1:.01:1.05;
    plot(x*ones(length(y)),y,'-k')
    hold on
end
ylim([0 1.2])
fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_all_areas_' mouse '_' date '_bouton_corr_strip.ps']);
        print(gcf, '-depsc2', fn_out);
        
%plot positions
ex_cells = [1 2 4 5 13 14];
ex_cells = [15 16 25 26 27 28];
ex_cells = [31 32 33 34 37 38];
ex_cells = [39 40 46 47 53 54];
ex_cells = [56 57 69 70 71 72];


col_ex = strvcat('ob','oc','or','om','og','oy'); 
figure;
dF_MAX = max(dF_nodir,[],3);
imagesq(dF_MAX);
colormap(gray)
caxis([0 .5])
hold on;
for iCell = 1:6
    scatter(pos(2,ex_cells(iCell)), pos(1,ex_cells(iCell)), 100, col_ex(iCell,:));
    hold on
end

%blank periods
fn_sorted =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted.tif']);
stack = readtiff(fn_sorted);

TC_sorted = [];
start= 1;
for ipair = 1:length(pairs)
    cell = order(:,pairs(ipair));
    iCell = goodfit_ind(:,cell);
    pos = all_fits(iArea).expt(iexp).bouton(iCell).pos;
    suby = pos(1)-1:pos(1)+1;
    subx = pos(2)-1:pos(2)+1;
    TC_sorted(:,start) = squeeze(mean(mean(stack(suby,subx,14305:end),2),1));
    start = start+1;
    if ipair == length(pairs)
        cell = order(:,pairs(ipair)+1);
        iCell = goodfit_ind(:,cell);
        pos = all_fits(iArea).expt(iexp).bouton(iCell).pos;
        suby = pos(1)-1:pos(1)+1;
        subx = pos(2)-1:pos(2)+1;
        TC_sorted(:,start) = squeeze(mean(mean(stack(suby,subx,14305:end),2),1));
        start = start+1;
    elseif pairs(ipair+1)>(pairs(ipair)+1)
        cell = order(:,pairs(ipair)+1);
        iCell = goodfit_ind(:,cell);
        pos = all_fits(iArea).expt(iexp).bouton(iCell).pos;
        suby = pos(1)-1:pos(1)+1;
        subx = pos(2)-1:pos(2)+1;
        TC_sorted(:,start) = squeeze(mean(mean(stack(suby,subx,14305:end),2),1));
        start = start+1;
    end   
end

TC_sorted_dF = bsxfun(@minus, TC_sorted,mean(TC_sorted,1));
TC_sorted_dFoverF = bsxfun(@rdivide, TC_sorted_dF, mean(TC_sorted,1));

figure; tcOffsetplot(TC_sorted_dFoverF(:,[1 2 4 5 13 14 15 16]));
figure; tcOffsetplot(TC_sorted_dFoverF(:,[25 26 27 28 31 32 33 34]));
figure; tcOffsetplot(TC_sorted_dFoverF(:,[37 38 39 40 46 47 53 54]));
figure; tcOffsetplot(TC_sorted_dFoverF(:,[56 57 69 70 71 72 73]));
hold on
for itime = 1:length(stim_times)
    x = 12:24:960;
    y = -2:1:20;
    plot(x(itime)*ones(length(y)),y,':k');
    hold on
end
ylim([-2 18]) 
fn_out = fullfile(fig_base, ['Figure' fig], [matrix '_' num2str(P) 'P_' inj '_all_areas_' mouse '_' date '_boutonTCs_blanks_row3.ps']);
        print(gcf, '-depsc2', fn_out);

