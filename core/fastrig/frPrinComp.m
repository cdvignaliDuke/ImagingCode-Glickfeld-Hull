function pcs = frPrinComp(newdir,varargin)
% FRPRINCOMP
% FRPRINCOMP(NEWDIR,VARARGIN)
% PCS = FRPRINCOMP(NEWDIR,...)
global stack
if ~sum(newdir==filesep)
    exptnames=frGetAllExpts(newdir)
    for index = 1:length(exptnames)
        frPrinComp(exptnames{index},varargin{:});
    end
    return;
end

[animal,exptname]=fileparts(newdir);

%%
defaultopts = {'nComp',300,'BorderWidth',4};
options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);

for iarg = 1:2:length(varargin)
    options.(varargin{iarg}) =varargin{iarg+1};
end;

%%
expt = frGetExpt(newdir);

% load stack
stack = single(readtiff(expt.dirs.reggreenpn));

% mask edges 
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
save(fullfile(expt.dirs.analrootpn,expt.filenames.pca_usv),'-struct','pcs');

return;
