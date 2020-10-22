function fast_errbar(xvals,data,trialDim,varargin)

p = inputParser;
p.addParamValue('shaded', false, @islogical);
p.addParamValue('color', [0 0 0], @isnumeric);

% parse inputs
p.parse(varargin{:});
params = p.Results;
    
mean_val = nanmean(data,trialDim);
std_val = nanstd(data,[],trialDim)/sqrt(size(data,trialDim));

if params.shaded
    h = shadedErrorBar(xvals,mean_val,std_val);
    h.mainLine.Color = params.color;
    h.patch.FaceColor = params.color;
else
    h = errorbar(xvals,nanmean(data,trialDim),nanstd(data,[],trialDim)/sqrt(size(data,trialDim)),'.','LineStyle','-','linewidth', 1.25,'MarkerSize',20);
    if ~isempty(params.color)
        h.Color = params.color;
    end
end