function ind1 = indgap(flag,varargin)
% function ind1 = indgap(flag,varargin) find the indexes corresponding of
% the start and end of each gap based on the flag assigned to the
% datapoints.

% Changes from the last version: gaps have flag==1
% Version: 1.0.3
% Author(s): Javier Pascual-Granado
% Date: 27/10/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(find(flag==1, 1)),
    ind1=[];
    return;
end

L = length(flag);
ind1 = nan(1,L);

% if nargin==1,
%     h = waitbar(0,'Step 0 - Finding indexes...');
% else
%     h = waitbar(0,sprintf('Step %d - Rebuilding indexes...',varargin{1}));
% end

% Main loop. Gap start and end indexes are saved in ind1
i=1;j=1;

while i <= L,
    % First index inside the gap
    ind = find(flag(i:end)==1,1);
    if ~isempty(ind),
        ind1(j) = i+ind-1;
        i = i+ind-1;
        j = j+1;
    else
        break;
    end
    
    % Last index inside the gap
    ind = find(flag(i:end)~= 1, 1) - 1;
    if ~isempty(ind),
        ind1(j) = i+ind-1;
        i = i+ind;
        j = j+1;
    else
        if flag(end)==1,
            ind1(j) = L;
        end
        break;
    end
%     waitbar(i/L,h);
end

if ~exist('ind1','var'),
    ind1=[];
    return;
end

% close(h);

% This sentence avoids some issues
ind1(isnan(ind1))=[];

end