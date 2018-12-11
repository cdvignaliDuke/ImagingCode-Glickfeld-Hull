% Size tuning fitting by KM, modified from Fit_2Dellipse_LG_Ret_KM

%params VB sends in: sta.kmean
% start with dumdum (dF/F responses for single cell+contrast)
% and szs0 (stim size for each response)

s = struct;
s.data = dumdum;
s.szs0 = szs0;

m1 = @(coefs,xdata) coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3))));
m2 = @(coefs,xdata) coefs(1)./(1+exp(-(coefs(2)+coefs(5))*(xdata-coefs(3)))) - coefs(4)./(1+exp(-coefs(5)*(xdata-(coefs(3)+coefs(6)))));
s.m1 = func2str(m1);
s.m2 = func2str(m2);

% ? what is CM
CM = zeros(1,2);
CM_data = zeros(1,2);
CM(1) = sum(x(:,1).*s.orig)/sum(s.orig);
CM(2) = sum(x(:,2).*s.orig)/sum(s.orig);
CM_data = CM;

% also not sure if mask_NaN needed, because all sizes have data points
Mask_NaN = isnan(y);
ind_NaN = find(isnan(y));
ind_noNaN = find(Mask_NaN == 0);
if ~isempty(ind_NaN)
    x(ind_NaN,:) = [];
    y(ind_NaN) = [];
end


% lower/upper bounds
%Ae ke xe
s.lb1 =  [0 0 0];
s.ub1 =  [inf 1 80];
%Ae k1 x1 Ai k2 x2
s.lb2 =  [0 0 0 0 0 0];
s.ub2 =  [inf 1 80 inf 1 80];

% for sampling different initial guesses
Nsamps = 5; %changed from 2 to 3
dbin1 = [s.ub1(1) - s.lb1(1); ...
    s.ub1(2) - s.lb1(2); ...
    s.ub1(3) - s.lb1(3)]./ (Nsamps-1);
dbin2 = [s.ub2(1) - s.lb2(1); ...
    s.ub2(2) - s.lb2(2); ...
    s.ub2(3) - s.lb2(3); ...
    s.ub2(4) - s.lb2(4); ...
    s.ub2(5) - s.lb2(5); ...
    s.ub2(6) - s.lb2(6)] ./ (Nsamps-1);

ke_vec = s.lb1(2)+dbin1(2)/2:dbin1(2):s.ub1(2);
xe_vec = s.lb1(3)+dbin1(3)/2:dbin1(3):s.ub1(3);

k1_vec = s.lb2(2)+dbin2(2)/2:dbin2(2):s.ub2(2);
x1_vec = s.lb2(3)+dbin2(3)/2:dbin2(3):s.ub2(3);
k2_vec = s.lb2(5)+dbin2(5)/2:dbin2(5):s.ub2(5);
x2_vec = s.lb2(6)+dbin2(6)/2:dbin2(6):s.ub2(6);

% loop over initial parameters not to get stuck in local minimum
% start with model1
clear temp;
index = 1;
names = {'c1','OF1','x0'};
% model1
for ike = 1:length(ke_vec)
    for ixe = 1:length(xe_vec)
        % fprintf('.');
        s.guess1 =  [0.5*maxMean ke_vec(ike) xe_vec(ixe)];
        
        [c1,OF1] = fminsearchbnd(@(c)sigevalpen1(c,s.szs0,s.data),s.guess1,s.lb1,s.ub1,opts);
        temp(index) = cell2struct({c1,OF1,s.guess1},names,2);
        index = index +1;
    end
end
[val,ind]=min([temp.OF1]);
s.fit1 = temp(ind);

% model2
clear temp;
index = 1;
names = {'c2','OF2','x0'};
for ik1 = 1:length(k1_vec)
    for ix1 = 1:length(x1_vec)
        for ik2 = 1:length(k2_vec)
            for ix2 = 1:length(x2_vec)
                % fprintf('.');
                s.guess2 =  [maxMean k1_vec(ik1) x1_vec(ix1) maxMean-mean(s.data) k2_vec(ik2) x2_vec(ix2)];
                
                [c2,OF2] = fminsearchbnd(@(c)sigevalpen2(c,s.szs0,s.data),s.guess2,s.lb2,s.ub2,opts);
                temp(index) = cell2struct({c2,OF2,s.guess2},names,2);
                index = index +1;
            end
        end
    end
end
[val,ind]=min([temp.OF2]);
s.fit2 = temp(ind);


% with selected fits, compute statistics
s.fitraw1 = m1(s.fit1.c1,s.szs0);
s.fitraw2 = m2(s.fit2.c2,s.szs0);
s.res1 = s.data-s.fitraw1;
s.res2 = s.data-s.fitraw2;
s.SStot = sum((s.data-mean(s.data)).^2);
s.SSE1 = norm(s.res1).^2; % SSE is squared 2-norm of error
s.SSE2 = norm(s.res2).^2;
s.Rsq1 = 1 - s.SSE1/s.SStot;
s.Rsq2 = 1 - s.SSE2/s.SStot;

szRng = linspace(0,max(szs)); % oversample at default 100 points
% find prefSize/SI for model1 (pref at 90% max)
s.fitout1 = m1(s.fit1.c1,szRng);
s.maxResp1 = max(s.fitout1);
s.prefSize1 = szRng(find(s.fitout1>(0.9*s.maxResp1),1));
if ~length(prefSize1)
    s.prefSize1=max(szs);
end
s.suppInd1 = 0;

% find prefSize/SI for model2 (pref at max)
s.fitout2 = m2(s.fit2.c2,szRng);
s.maxResp2 = max(s.fitout2);
s.prefSize2 = szRng(find(s.fitout2==s.maxResp2,1));
s.suppInd2 = 1 - s.fitout2(end)/s.maxResp2;


% plot both fits (with initial guesses) and residual (4 plots each, in blocks)
if PLOTIT_FIT == 1
    if start ==64
        set(gcf, 'Position', [0 0 800 1000]);
        fn_out = fullfile('\\CRASH.dhe.duke.edu\data\home\kevin\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_RFfits' num2str(ifig) '.pdf']);
        print(fn_out,'-dpdf')
        figure;
        ifig = 1+ifig;
        start = 1;
    end
    h = subplot(8,8,start);
    errorbar([0 szs],[0 sizeMean(:,nCon,iCell)'],[0 sizeSEM(:,nCon,iCell)'])
    hold on
    plot(s.szs0,s.data,'.b')
    plot(szRng,s.fitout1,'-')
    plot(szRng,m1(s.fit1.x0,szRng),'g--')
    hold off
    ylim([min([-0.5*s.maxResp1 min(sizeMean(:,nCon,iCell))]) 1.4*max([s.maxResp1 max(sizeMean(:,nCon,iCell))])])
    title([num2str(iCell), ' prefsz1 ' num2str(s.prefSize1)]);
    
    h = subplot(8,8, 1+start);
    plot(s.szs0,s.res1)
    title(['R^2_1: ' num2str(s.Rsq1)])
    
    h = subplot(8,8, 2+start);
    errorbar([0 szs],[0 sizeMean(:,nCon,iCell)'],[0 sizeSEM(:,nCon,iCell)'])
    hold on
    plot(s.szs0,s.data,'.b')
    plot(szRng,s.fitout2,'-')
    plot(szRng,m2(s.fit2.x0,szRng),'g--')
    hold off
    ylim([min([-0.5*s.maxResp2 min(sizeMean(:,nCon,iCell))]) 1.4*max([s.maxResp2 max(sizeMean(:,nCon,iCell))])])
    title([num2str(iCell), ' prefsz2 ' num2str(s.prefSize2)]);
    
    h = subplot(5,5,3+start);
    plot(s.szs0,s.res2)
    title(['R^2_2: ' num2str(s.Rsq2)])
    
    start = start+4;
end

if SAVEALLDATA == 0
    s = rmfield(s,'data');
    s = rmfield(s,'guess1');
    s = rmfield(s,'guess2');
    
end
