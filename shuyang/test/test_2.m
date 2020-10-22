%% outputs:
%   c: T x 1 vector, denoised trace
%   s: T x 1 vector, deconvolved signal
%   b: fluorescence baseline
%   kernel: struct variable containing the parameters for the selected
%       convolution model
%   lambda: Optimal Lagrange multiplier for noise constraint under L1 penalty
%     """olves the noise constrained sparse nonnegat

%%
c = [];
s = [];
threshold_stay = -4;
for i = 12
    [c, s, options] = deconvolveCa(TCave(:,i), 'optimize_pars', true, ...
        'optimize_b', true, 'method','foopsi', 'smin', threshold_stay);
    %[kernel_stay(:,c), spk_stay(:,c), options] = deconvolveCa();
end

figure;
subplot(3,1,1);
plot(TCave(:,12));
ylabel('rawF');
subplot(3,1,2);
plot(c); hold on; plot(s);
ylabel('denoised signal overlay with deconvolved signal');


subplot(3,1,3);
plot(s);
ylabel('deconvolved signal');
xlabel('frames');
