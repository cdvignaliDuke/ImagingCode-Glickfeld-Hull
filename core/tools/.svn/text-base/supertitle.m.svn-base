function f = supertitle(s,pos,varargin)
% SUPERTITLE 	makes a big title over all subplots
%
%	supertitle(S) writes the string S as a title
%	It returns the handle if you want it.
%
%	supertitle(S,pos) lets you specify a position be 0 and 1. Default is .95
%
% 09/07/01 vb re-implementat with annotation (to preserve axis ordering)
% 10/03/27 vb changed default pos to 0.99. 
%             varargin now passed to annotation
%             default backgroundcolor is w

if nargin < 2 
    pos = [];
end
    
if isempty(pos);pos = .99; end 

f = annotation('textbox',[0 pos 1 .01],'string',s,'edgecolor','none','hori','center',...
               'background','w',varargin{:},'interpreter','none');
return;

% currax = gca;
% 
% h = axes('position',[.5 pos 0.1 0.1 ],'xcolor',[0 0 0],'visible','off');
% 
% set(h,'color','none');
% 
% if ~isempty(s)
% 	f = text(0 , 0, s, 'units', 'normalized',...
% 		'verticalalignment','top','horizontalalignment','center' );
% end
% 
% axes(currax);
