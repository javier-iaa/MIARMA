function seg = maxseg(s, stat, varargin)
% function seg = maxseg(s, stat, varargin)
% This gives the largest segment without gaps. Inputs are the time series s 
% and the corresponding status array stat which is 1 for gaps (zero-valued 
% in s) and 0 otherwise.
%
% The output is the segment of max length. Also if the function is called
%     maxseg(s, stat, 'w', fname)
% the output is written to a file named by the string fname.
%
% Note: use regsamp.m to obtain the input s.
% 
% By Javier Pascual-Granado
% <a href="matlab:web http://www.iaa.csic.es;">IAA-CSIC, Spain</a>
%
% Version: 0.1
%
% Date: 14/04/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    flag = stat;

    if isempty( find(flag~=0, 1) )
        seg = s;
        return
    end

    L = length(flag);
    igap = nan(1,L);

    i=1;
    j=1;

    % Determine gap indexes
    while i <= L,
        ind = find(flag(i:end) ==1,1);
        if ~isempty(ind),
            igap(j) = i+ind-1;
            i = i+ind-1;
            j = j+1;
        else
            break;
        end

        ind = find(flag(i:end) == 0,1)-1;
        if ~isempty(ind),
            igap(j) = i+ind-1;
            i = i+ind;
            j = j+1;
        else
            break;
        end
    end

    segl = [igap(1)-1 (igap(3:2:end-1)-igap(2:2:end-2)-1) L-igap(end)];
    ML = max(segl);

    % The index of the corresponding segment
    I = find(segl==ML,1);

    % The largest segment
    if I==1
        seg = s(1:(igap(1)-1));
    elseif I==length(segl)
        seg = s( (igap(end)+1):end ) ;
    else
        seg = s((igap(I*2-2)+1):(igap(I*2-1)-1));
    end
    
    if nargin>2
        if strcmp(varargin{1}, 'w')
            fname = varargin{2};
            fich = fopen( fname, 'w');
            for i=1:length(seg)
                fprintf(fich, '%16.13f\n', seg(i));
            end
            fclose(fich);
        end
    end

end
