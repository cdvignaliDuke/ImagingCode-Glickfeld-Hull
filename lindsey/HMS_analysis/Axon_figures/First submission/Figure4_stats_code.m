fig_base = '\\zmey\storlab\users\Lindsey\Projects\HVAs\Manuscript\Figures_2012';
fig = 4;
col = strvcat('c', 'k', 'r');
P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';

fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
load(fn_good);

vec0 = zeros(2,7);
vec0(1,:) = [.5 1 2 4 8 16 32];
vec0(2,:) = [.01 .02 .04 .08 .16 .32 .64];    
vec1 = interp2(vec0');
TF_vec = zeros(8,1);
SF_vec = zeros(8,1);
TF_vec(1,1)= vec1(2,1);
TF_vec(2,1)= 1.0001;
TF_vec(3:6,1) = vec1(4:2:10,1);
TF_vec(7,1) = 14.999;
TF_vec(8,1) = vec1(12,1);
SF_vec(1,1)= vec1(2,3);
SF_vec(2,1)= .020001;
SF_vec(3:6,1) = vec1(4:2:10,3);
SF_vec(7,1) = .3199;
SF_vec(8,1) = vec1(12,3);
edges = log2(TF_vec./flipud(SF_vec));
plot_edges = [edges(2) mean(edges(2:3)) mean(edges(3:4)) mean(edges(4:5)) mean(edges(5:6)) mean(edges(6:7)) edges(7)];
SF_plot_edges = log2([.02 .025 .045 .09 .18 .28 .32]);
TF_plot_edges = log2([1 1.25 2.25 4.5 9 13.5 15]);

%comparison of speeds less than 25 and more than 100 deg/sec
p_dF = zeros(1,3);
for iArea = 1:3
    [n_speed bin_speed] = histc(log2(Goodfits(iArea).speed), edges);
    slow_speed_ind = find(bin_speed <4);
    fast_speed_ind = find(bin_speed >4);
    slow_dF = Goodfits(iArea).dF(slow_speed_ind);
    fast_dF = Goodfits(iArea).dF(fast_speed_ind);
    [h_dF, p_dF(:,iArea)] = ttest2(slow_dF, fast_dF);
end

%comparison of speeds in AL and PM in different dF/F bins
p_AL_PM_dF_speed = zeros(1,5);
[n_dF_PM bin_dF_PM] = histc(Goodfits(1).dF, edges_dF);
[n_dF_AL bin_dF_AL] = histc(Goodfits(3).dF, edges_dF);
for ibin = 1:length(n_dF)-1
    ind_PM = find(bin_dF_PM == ibin);
    ind_AL = find(bin_dF_AL == ibin);
    [h_AL_PM_dF_speed, p_AL_PM_dF_speed(1,ibin)] = kstest2(log2(Goodfits(1).speed(ind_PM)), log2(Goodfits(3).speed(ind_AL)));
end

%discriminability index
n_speed_AL_PM = zeros(8,2,5);
[n_dF_PM bin_dF_PM] = histc(Goodfits(1).dF, edges_dF);
[n_dF_AL bin_dF_AL] = histc(Goodfits(3).dF, edges_dF);
for ibin = 1:length(n_dF)-1
    ind_PM = find(bin_dF_PM == ibin);
    ind_AL = find(bin_dF_AL == ibin);
    [n_speed_PM bin_speed_PM] = histc(log2(Goodfits(1).speed(ind_PM)), edges);
    [n_speed_AL bin_speed_AL] = histc(log2(Goodfits(3).speed(ind_AL)), edges);
    n_speed_AL_PM(:,:,ibin) = [n_speed_AL n_speed_PM];
end

ErrMin_TF = zeros(1,5);
ErrMin_SF = zeros(1,5);
ErrMin_speed = zeros(1,5);
for ibin = 1:length(n_dF)-1
    ind_PM = find(bin_dF_PM == ibin);
    ind_AL = find(bin_dF_AL == ibin);

    SF_vec1 = log2(Goodfits(3).SF(ind_AL));
    SF_vec2 = log2(Goodfits(1).SF(ind_PM));
    TF_vec1 = log2(Goodfits(3).TF(ind_AL));
    TF_vec2 = log2(Goodfits(1).TF(ind_PM));
    speed_vec1 = log2(Goodfits(3).speed(ind_AL));
    speed_vec2 = log2(Goodfits(1).speed(ind_PM));

    SL = [SF_vec1; SF_vec2];
    SW = [TF_vec1; TF_vec2];
    SX = [speed_vec1; speed_vec2];

    group = [1*ones(length(SF_vec1),1); 2*ones(length(SF_vec2),1)];

    TF2 = SW;
    [i,j] = sort(TF2);
    val1 = (i(1:(end-1)) + i(2:end))/2;
    group1 = group(j);
    PercErr = zeros(size(val1));
    for count1 = 1:length(PercErr)
       val11 = val1(count1);
       ind1 = find(i<val11);
       N1_correct = length(find(i<val11 & group1==2));
       N1_false = length(find(i<val11 & group1==1));
       N2_correct = length(find(i>=val11 & group1==1));
       N2_false = length(find(i>=val11 & group1==2));
       PercErr(count1) = (N1_false + N2_false)./(length(i));
    end
    [i2,j2] = min(PercErr);
    ErrMin_TF(1,ibin) = 1-i2;
    ErrMin_TF_val = val1(j2);
    SF2 = SL;
    [i,j] = sort(SF2);
    val1 = (i(1:(end-1)) + i(2:end))/2;
    group1 = group(j);
    PercErr = zeros(size(val1));
    for count1 = 1:length(PercErr)
       val11 = val1(count1);
       ind1 = find(i<val11);
       N1_correct = length(find(i<val11 & group1==1));
       N1_false = length(find(i<val11 & group1==2));
       N2_correct = length(find(i>=val11 & group1==2));
       N2_false = length(find(i>=val11 & group1==1));
       PercErr(count1) = (N1_false + N2_false)./(length(i));
    end
    [i2,j2] = min(PercErr);
    ErrMin_SF(1,ibin) = 1-i2;
    ErrMin_SF_val = val1(j2);
    speed2 = SX;
    [i,j] = sort(speed2);
    val1 = (i(1:(end-1)) + i(2:end))/2;
    group1 = group(j);
    PercErr = zeros(size(val1));
    for count1 = 1:length(PercErr)
       val11 = val1(count1);
       ind1 = find(i<val11);
       N1_correct = length(find(i<val11 & group1==2));
       N1_false = length(find(i<val11 & group1==1));
       N2_correct = length(find(i>=val11 & group1==1));
       N2_false = length(find(i>=val11 & group1==2));
       PercErr(count1) = (N1_false + N2_false)./(length(i));
    end
    [i2,j2] = min(PercErr);
    ErrMin_speed(1,ibin) = 1-i2;
    ErrMin_speed_val = val1(j2);
end


%comparisons of SF/TF/speed for boutons with dF less than 0.4 and more than 0.8
[n_dF_PM bin_dF_PM] = histc(Goodfits(1).dF, edges_dF);
[n_dF_AL bin_dF_AL] = histc(Goodfits(3).dF, edges_dF);
low_dF_PM_ind = find(bin_dF_PM<3);
high_dF_PM_ind = find(bin_dF_PM>3);
low_dF_AL_ind = find(bin_dF_AL<3);
high_dF_AL_ind = find(bin_dF_AL>3);

[h_speed_PM, p_speed_PM] = ttest2(log2(Goodfits(1).speed(low_dF_PM_ind)), log2(Goodfits(1).speed(high_dF_PM_ind)));
[h_SF_PM, p_SF_PM] = ttest2(log2(Goodfits(1).SF(low_dF_PM_ind)), log2(Goodfits(1).SF(high_dF_PM_ind)));
[h_TF_PM, p_TF_PM] = ttest2(log2(Goodfits(1).TF(low_dF_PM_ind)), log2(Goodfits(1).TF(high_dF_PM_ind)));

[h_speed_AL, p_speed_AL] = ttest2(log2(Goodfits(3).speed(low_dF_AL_ind)), log2(Goodfits(3).speed(high_dF_AL_ind)));
[h_SF_AL, p_SF_AL] = ttest2(log2(Goodfits(3).SF(low_dF_AL_ind)), log2(Goodfits(3).SF(high_dF_AL_ind)));
[h_TF_AL, p_TF_AL] = ttest2(log2(Goodfits(3).TF(low_dF_AL_ind)), log2(Goodfits(3).TF(high_dF_AL_ind)));


    
    