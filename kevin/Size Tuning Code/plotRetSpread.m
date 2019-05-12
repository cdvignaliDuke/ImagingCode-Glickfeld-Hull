% receptive field spread
clear all; clc;
%expfile = '\\CRASH.dhe.duke.edu\data\home\kevin\Code\Ai9x_experiment_list.txt';
expfile = 'C:\Users\kevin\Documents\Repositories\ImagingCode-Glickfeld-Hull\kevin\Size Tuning Code\Ai9x_experiment_list.txt';
fID = fopen(expfile);
head = textscan(fID,'%s%s%s%s%s',1,'delimiter',',');
head = vertcat(head{:});
temp = textscan(fID,'%s%s%s%s%s','delimiter',',','HeaderLines',1);
temp = horzcat(temp{:});
expdata = cell2table(temp,'VariableNames',head);
nExp = size(expdata,1);
%isvalid = ones(1,nExp);
%expdata = addvars(expdata,isvalid);

fprintf(['Receptive-field spread figure - by KM, Glickfeld Lab\nLoading ' num2str(nExp) ' experiments\n'])

%% load each experiment and concatenate data

% load ret exp data
fprintf('Loading lbub_fits\n')
lbub_fits_all = [];
goodfit_ind_all = [];
cellAz_norm_all = [];
cellEl_norm_all = [];
nCellsExpRet = zeros(1,nExp);
for i=1:nExp
    fprintf(['Exp: ' num2str(i) '/' num2str(nExp) '...'])
    
    date = expdata.date{i};
    mouse = expdata.mouse{i};
    run_str = 'runs-002'; %expdata.run_str{i};
    load(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOVs.mat']), 'Azs', 'Els')
    filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
    if ~exist(filename, 'file')
        fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
    end
    load(filename, 'lbub_fits', 'goodfit_ind')
    nCellsExpRet(i) = size(lbub_fits,1);
    
    if i==49%[25 39])
        figure(1);clf;
        for i2 = goodfit_ind
            plot(lbub_fits(i2,4,4), lbub_fits(i2,5,4), 'ok')
            hold on;
            if i2==47
                text(lbub_fits(i2,4,4), lbub_fits(i2,5,4), num2str(i2),'Color','r','HorizontalAlignment','center','VerticalAlignment','middle')
            end
        end
        xlim([min(Azs,[],2) max(Azs,[],2)])
        ylim([min(Els,[],2) max(Els,[],2)]+5)
        axis square
        %title('RF center')
        xlabel('Azimuth (deg)')
        ylabel('Elevation (deg)')
        
        figure(2);clf;
        for i2 = goodfit_ind
            if i2==47
                h = ellipse(sqrt(2*log(2))*lbub_fits(i2,2,4), sqrt(2*log(2))*lbub_fits(i2,3,4), 0, lbub_fits(i2,4,4), lbub_fits(i2,5,4),'r');
                %h.Color = [1 0 0 0.1];
            else
                h = ellipse(sqrt(2*log(2))*lbub_fits(i2,2,4), sqrt(2*log(2))*lbub_fits(i2,3,4), 0, lbub_fits(i2,4,4), lbub_fits(i2,5,4));
                h.Color = [0 0 0 0.1];
            end
            hold on;
        end
        xlim([min(Azs,[],2) max(Azs,[],2)])
        ylim([min(Els,[],2) max(Els,[],2)]+5)
        axis square
        %title('RF center')
        xlabel('Azimuth (deg)')
        ylabel('Elevation (deg)')
        pause
    end
    
    fprintf('done\n')
end

fprintf(['\nFinished loading all ' num2str(nExp) ' experiments.\n'])

%% load a selected experiment
iExp = 49;
fprintf(['Loading exp: ' num2str(iExp) '/' num2str(nExp) '...'])

date = expdata.date{iExp};
mouse = expdata.mouse{iExp};
run_str = 'runs-002'; %expdata.run_str{i};
load(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_FOVs.mat']), 'Azs', 'Els')
filename = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_lbub_fits.mat']);
if ~exist(filename, 'file')
    fprintf([[date '_' mouse '_' run_str '_lbub_fits.mat'] ' not found! Please remove from list\n'])
end
load(filename, 'lbub_fits', 'goodfit_ind')
nCellsExpRet(iExp) = size(lbub_fits,1);

figure(1);clf;
for i2 = goodfit_ind
    if i2==47
        %text(lbub_fits(i2,4,4), lbub_fits(i2,5,4), num2str(i2),'Color','r','HorizontalAlignment','center','VerticalAlignment','middle')
        plot(lbub_fits(i2,4,4), lbub_fits(i2,5,4), 'or', 'MarkerFaceColor', 'r')
    elseif i2==254
        plot(lbub_fits(i2,4,4), lbub_fits(i2,5,4), 'ob', 'MarkerFaceColor', 'b')
    else
        plot(lbub_fits(i2,4,4), lbub_fits(i2,5,4), '.k')
    end
    hold on
end
xlim([min(Azs,[],2) max(Azs,[],2)])
ylim([min(Els,[],2) max(Els,[],2)])
axis square
%title('RF center')
xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')

figure(2);clf;
for i2 = goodfit_ind
    if i2==47
        h = ellipse(sqrt(2*log(2))*lbub_fits(i2,2,4), sqrt(2*log(2))*lbub_fits(i2,3,4), 0, lbub_fits(i2,4,4), lbub_fits(i2,5,4),'r');
        %h.Color = [1 0 0 0.1];
    elseif i2==254
        h = ellipse(sqrt(2*log(2))*lbub_fits(i2,2,4), sqrt(2*log(2))*lbub_fits(i2,3,4), 0, lbub_fits(i2,4,4), lbub_fits(i2,5,4),'b');
    else
        h = ellipse(sqrt(2*log(2))*lbub_fits(i2,2,4), sqrt(2*log(2))*lbub_fits(i2,3,4), 0, lbub_fits(i2,4,4), lbub_fits(i2,5,4));
        h.Color = [0 0 0 0.1];
    end
    hold on;
end
xlim([min(Azs,[],2) max(Azs,[],2)])
ylim([min(Els,[],2) max(Els,[],2)])
axis square
set(gca,'box','on','TickDir','in')
%title('RF center')
xlabel('Azimuth (deg)')
ylabel('Elevation (deg)')

fprintf('done\n')

%% present two cells
% load _Tuning.mat: 'tc_dfof', 'tuning_mat', 'Stims', 'Ind_struct'
load(fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Tuning.mat']))


Azs00 = min(Azs):0.1:max(Azs);
Els00 = max(Els):-0.1:min(Els);
[p q] = meshgrid(Azs00,Els00); Stims00 = [p(:) q(:)];

figure(3);clf;
set(gcf, 'Position', [100 100 520 510]); %[100 100 650 640]

% first cell (47)
iCell=47;
subplot(2,2,1)
ret_raw = reshape(tuning_mat(:,1,iCell), [length(Azs) length(Els)])';
dFoFmax = max(ret_raw(:));
ret_norm = (ret_raw+min(ret_raw(:)))./(dFoFmax+min(ret_raw(:)));
ret_red = cat(3,ret_norm,0*ret_norm,0*ret_norm);
imagesc(ret_red)
title(['Max dF/F = ' num2str(dFoFmax,3)])
axis equal
xticks(1:7); xticklabels(Azs); yticks(1:7); yticklabels(Els)
xlabel('Azimuth (deg)'); ylabel('Elevation (deg)')
% fit
subplot(2,2,2)
fit_mat = Gauss2D_ellipseMA(lbub_fits(iCell,:,4),Stims00);%get gaussian as 25x1 vector
ret_fit = reshape(fit_mat, [length(Azs00) length(Els00)]);
dFoFmax = max(ret_fit(:));
ret_norm = (ret_fit+min(ret_fit(:)))./(dFoFmax+min(ret_fit(:)));
ret_red = cat(3,ret_norm,0*ret_norm,0*ret_norm);
imagesc(Azs00,Els00,ret_red)
set(gca,'YDir','normal')
hold on
%plot(lbub_fits(iCell,4,4), lbub_fits(iCell,5,4),'wx')
errorbar(lbub_fits(iCell,4,4), lbub_fits(iCell,5,4), lbub_fits(iCell,5,4)-lbub_fits(iCell,5,1), lbub_fits(iCell,5,2)-lbub_fits(iCell,5,4), lbub_fits(iCell,4,4)-lbub_fits(iCell,4,1), lbub_fits(iCell,4,2)-lbub_fits(iCell,4,4),'w')
ellipse(sqrt(2*log(2))*lbub_fits(iCell,2,4), sqrt(2*log(2))*lbub_fits(iCell,3,4), 0, lbub_fits(iCell,4,4), lbub_fits(iCell,5,4),'w');
colormap gray
title(['2D Gaussian fit'])
axis equal
xlabel('Azimuth (deg)'); ylabel('Elevation (deg)')
xlim([min(Azs,[],2) max(Azs,[],2)]); ylim([min(Els,[],2) max(Els,[],2)])

% second cell (254)
iCell = 254;
subplot(2,2,3)
ret_raw = reshape(tuning_mat(:,1,iCell), [length(Azs) length(Els)])';
dFoFmax = max(ret_raw(:));
ret_norm = (ret_raw+min(ret_raw(:)))./(dFoFmax+min(ret_raw(:)));
ret_blue = cat(3,0*ret_norm,0*ret_norm,ret_norm);
imagesc(ret_blue)
colormap gray
title(['Max dF/F = ' num2str(dFoFmax,2)])
axis equal
xticks(1:7); xticklabels(Azs); yticks(1:7); yticklabels(Els)
xlabel('Azimuth (deg)'); ylabel('Elevation (deg)')
% fit
subplot(2,2,4);cla;
fit_mat = Gauss2D_ellipseMA(lbub_fits(iCell,:,4),Stims00);%get gaussian as 25x1 vector
ret_fit = reshape(fit_mat, [length(Azs00) length(Els00)]);
ret_fit = ret_fit;
dFoFmax = max(ret_fit(:));
ret_norm = (ret_fit+min(ret_fit(:)))./(dFoFmax+min(ret_fit(:)));
ret_blue = cat(3,0*ret_norm,0*ret_norm,ret_norm);
imagesc(Azs00,Els00,ret_blue)
set(gca,'YDir','normal')
hold on
%plot(lbub_fits(iCell,4,4), lbub_fits(iCell,5,4),'wx')
errorbar(lbub_fits(iCell,4,4), lbub_fits(iCell,5,4), lbub_fits(iCell,5,4)-lbub_fits(iCell,5,1), lbub_fits(iCell,5,2)-lbub_fits(iCell,5,4), lbub_fits(iCell,4,4)-lbub_fits(iCell,4,1), lbub_fits(iCell,4,2)-lbub_fits(iCell,4,4),'w')
ellipse(sqrt(2*log(2))*lbub_fits(iCell,2,4), sqrt(2*log(2))*lbub_fits(iCell,3,4), 0, lbub_fits(iCell,4,4), lbub_fits(iCell,5,4),'w');
colormap gray
title(['2D Gaussian fit'])
axis equal
xlabel('Azimuth (deg)'); ylabel('Elevation (deg)')
xlim([min(Azs,[],2) max(Azs,[],2)]); ylim([min(Els,[],2) max(Els,[],2)])

%% try to do timecourses?
figure(1);clf;
iCell1=47; %47, 254
iCell2 = 254;
h = zeros(49,1);
max_tc1=0; max_tc2=0;
tc_cell1 = zeros(size(tc_dfof,1),length(Stims));
tc_cell2 = tc_cell1;
time = [1:15]*1/15*1000
for i=1:length(Stims)
    tc_cell1(:,i) = mean(tc_dfof(:,iCell1,[Ind_struct(i).all_trials]),3);
    tc_cell2(:,i) = mean(tc_dfof(:,iCell2,[Ind_struct(i).all_trials]),3);
    max_tc1 = max(max_tc1,max(tc_cell1(:,i)));
    max_tc2 = max(max_tc2,max(tc_cell2(:,i)));
end
% normalize
tc_cell1 = tc_cell1./max_tc1;
tc_cell2 = tc_cell2./max_tc2;
for i=1:length(Stims)
    h(i)=subplot(7,7,i);
    plot(time,tc_cell1(46:end,i),'r');
    hold on
    plot(time,tc_cell2(46:end,i),'b');
    set(gca,'box','off','TickDir','out')
    if i==22; ylabel('dF/F (norm)'); end
    if i==46; xlabel('Time (ms)'); end
end
linkaxes(h,'xy')
max_tc = 1% max(max_tc1,max_tc2);
ylim([-0.2*max_tc 1.2*max_tc])