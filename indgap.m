function ind1 = indgap(flag)
% function ind1 = indgap(flag) find the indexes corresponding of
% the start and end of each gap based on the flag assigned to the
% datapoints.

% Changes from the last version: minor changes
% Version: 1.0.5
% Author(s): Javier Pascual-Granado
% Date: 28/04/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(find(flag==1, 1)),
    ind1=[];
    return;
end

L = length(flag);
ind1 = nan(1,L);

% Main loop. Gap start and end indexes are saved in ind1
i=1;
j=1;

while i <= L,
    % First index inside the gap
    ind = find( abs( flag(i:end) )==1, 1);
    if ~isempty(ind),
        ind1(j) = i+ind-1;
        i = i+ind-1;
        j = j+1;
    else
        break;
    end
    
    % Last index inside the gap
    ind = find(abs(flag(i:end))~= 1, 1) - 1;
    if ~isempty(ind),
        ind1(j) = i+ind-1;
        i = i+ind;
        j = j+1;
    else
        if abs(flag(end))==1,
            ind1(j) = L;
        end
        break;
    end

end

if ~exist('ind1','var'),
    ind1=[];
    return;
end

% This sentence avoids some issues
ind1(isnan(ind1))=[];

end
