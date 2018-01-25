function [y_fit, theta_dist_save,theta_90,R_square] = vonmisesReliableFit_direction(data, data_resamp,theta,nboot)

% data is the average response to each stimulus for each cell: [nCells, nOri+1] 0:180/nOri:180
% data_resamp the average response for each cell randomly resampling trials
% 1000 times [nCells, nOrii+1, nboot]
nCells = size(data,1);
theta_smooth = 0:1:360;
b = nan(nboot+1,nCells);
k = nan(nboot+1,nCells);
R1 = nan(nboot+1,nCells);
u1 = nan(nboot+1,nCells);
R2 = nan(nboot+1,nCells);
u2 = nan(nboot+1,nCells);
R_square = nan(nboot+1,nCells);
y_fit = nan(length(theta_smooth),nboot+1,nCells);
max_ori = nan(nboot+1,nCells);

for iCell = 1:nCells
    if all(data(iCell,:) < 0)
        y_fit(:,1,iCell) = zeros(size(y_fit,1),1,1);
        y_fit(:,2:nboot+1,iCell) = nan(size(y_fit,1),nboot,1);
    else        
        [b(1,iCell), k(1,iCell), R1(1,iCell), R2(1,iCell),...
            u1(1,iCell), u2(1,iCell),~,R_square(1,iCell)] = ...
            miaovonmisesfit_dir(deg2rad(theta),data(iCell,:));
        y_fit(:,1,iCell) = b(1,iCell)+R1(1,iCell).*exp(k(1,iCell)...
                .*(cos((deg2rad(theta_smooth)-u1(1,iCell)))-1)) + ...
                b(1,iCell)+R2(1,iCell).*exp(k(1,iCell)...
                .*(cos((deg2rad(theta_smooth)-u2(1,iCell)))-1));
        [~, max_ori(1,iCell)] = max(y_fit(:,1,iCell),[],1);
        for i = 2:nboot+1
            [b(i,iCell), k(i,iCell), R1(i,iCell), R2(i,iCell),...
            u1(i,iCell), u2(i,iCell),~,R_square(i,iCell)] = ...
            miaovonmisesfit_dir(deg2rad(theta),data_resamp(iCell,:,i-1));
            y_fit(:,i,iCell) = b(i,iCell)+R1(i,iCell).*exp(k(i,iCell)...
                .*(cos((deg2rad(theta_smooth)-u1(i,iCell)))-1)) + ...
                b(i,iCell)+R2(i,iCell).*exp(k(i,iCell)...
                .*(cos((deg2rad(theta_smooth)-u2(i,iCell)))-1));
            [~, max_ori(i,iCell)] = max(y_fit(:,i,iCell),[],1);
        end
        fprintf('.')
    end
end

theta_90 = nan(1,nCells);
theta_dist_save = nan(nboot,nCells);
for iCell = 1:nCells
    if all(data(iCell,:) < 0)
        theta_dist_save(1,iCell) = nan;
        theta_90(1,iCell) = nan;
    else
        if ~isnan(R_square(iCell))
            theta_dist = abs(theta_smooth(squeeze(max_ori(1,iCell)))...
                -theta_smooth(squeeze(max_ori(2:nboot+1,iCell))));
            theta_dist(find(theta_dist>90)) = 180-theta_dist(find(theta_dist>90));
            theta_dist_save(:,iCell) = theta_dist;
            theta_sort = sort(theta_dist,'ascend');
            theta_90(1,iCell) = theta_sort(ceil(nboot*.9));
        end
    end
end
end