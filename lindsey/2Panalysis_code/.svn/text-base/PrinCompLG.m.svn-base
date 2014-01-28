function pcs = PrinCompLG(fn)
% FRPRINCOMP
% FRPRINCOMP(NEWDIR,VARARGIN)
% PCS = FRPRINCOMP(NEWDIR,...)
global stack

defaultopts = {'nComp',300,'BorderWidth',4};
options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);


% mask edges 
stack = single(readtiff(fn));
[ny,nx,nt]=size(stack);
roi = imborder([ny,nx],options.BorderWidth,0); 
fprintf('Masking edges... ');
stack= bsxfun(@times,stack,single(roi));
fprintf('Done\n');

% compute thin svd using randomized algorithm
%pcs = stackFastPCA(stack,options.nComp);
pcs = stackFastPCA(1,options.nComp);
% save extracted components 
fprintf('Saving principal components\n');
fn_out = [fn(1:size(fn,2)-4) '_pca_usv.mat'];
save(fn_out,'-struct','pcs');

return;