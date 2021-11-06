function [data1, flag1] = lincorr(data, flag, ind, npi)
% function [data1,flag1] = lincorr(data,flag,ind,npi) use a polynomial
% interpolation for filling the smallest gaps in the time series.
% Inputs:   data - time series
%               flag - status array
%               ind - gap indexes
%               npi - limits the size of the gaps to be interpolated.
%
% Ouputs:  data1 - time series after lin.interpolation
%                flag1 - new status array
%
% Version: 1.1
%
% Changes from the last version:
% - Code completely refurbished.
% - Now using previous interpolated data points
%
% Calls: polintre.m
% Author(s): Javier Pascual-Granado
% Date: 17/05/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Order used for the polynomial interpolation
ord = 3;

% Data points used for the polynomial fitting. Be careful with this, just
% a few data points is recommended.
npint = ord*2;

if npi>npint
    npint = npi;
end

flag1 = flag;
data1 = data;

N = length(data);

% Change column into row
data1 = reshape(data1,1,N);

if isempty(ind)==1,
    return;
end

L = length(ind);
last_gap_ind = ind(end-1);

% Initialization
ini_gap_ind = ind(1);
end_gap_ind = ind(2);
if L>2
    next_gap_ind = ind(3);
end

while ini_gap_ind <= last_gap_ind
    
    dif = end_gap_ind - ini_gap_ind + 1;
    
    if dif <= npi
        % Select left segment
        subi1 = 1:ini_gap_ind - 1;
        is_gap = find( flag1(subi1)==1 );
        if ~isempty( is_gap )
            subi1 = subi1( is_gap(end)+1: end);
        end
        
        if length(subi1) > npint
            subi1(1:end - npint) = [];
        end
        
        seg1 = data1( subi1 );
        
        % Select right segment
        if (L==2) || (ini_gap_ind==last_gap_ind)
            subi2 = end_gap_ind + 1:N;
        else
            subi2 = end_gap_ind + 1:next_gap_ind - 1;
        end

        if length(subi2) > npint
            subi2 = subi2(1:npint);
        end
        
        seg2 = data1( subi2 );

        if isempty(seg1) || isempty(seg2),
            continue;
        end
   
        interp = polintre (seg1, seg2, dif, ord);
        subint = ini_gap_ind:end_gap_ind;
        data1( subint ) = interp;
        flag1( subint ) = 0;
    end

    if (L==2) || (ini_gap_ind==last_gap_ind)
        break;
    else    
        ind(1:2) = [];
        ini_gap_ind = next_gap_ind;
        end_gap_ind = ind(2);
        if next_gap_ind==last_gap_ind
            next_gap_ind = [];
        else
            next_gap_ind = ind(3);
        end
    end
end
