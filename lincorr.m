function [data1, flag1] = lincorr(data, flag, ind, npi)
% function [data1,flag1] = lincorr(data,flag,ind,npi) use a linear
% interpolation for gap-filling the smallest gaps. 
% Inputs:   data - time series
%           flag - status array
%           ind - gap indexes
%           npi - limits the size of the gaps to be interpolated.
% Ouputs:   data1 - time series after lin.interpolation
%           flag1 - new status array
% Version: 1.0.7
% Changes from the last version:
%   - Code cleaned.
%   - Calls polintre v0.1
%   - Fixed issue with npint being too small.
%
% Calls: polintre.m
% Author(s): Javier Pascual-Granado
% Date: 05/03/2021
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

i=1;
while i <= L-1,
    if L==2,
        seg1 = data1(1:(ind(1)-1));
        seg2 = data1((ind(2)+1):end);
    elseif i==1,
        seg1 = data1(1:(ind(i)-1));
        seg2 = data1((ind(i+1)+1):(ind(i+2)-1));
    elseif i==L-1,
        seg1 = data1((ind(L-2)+1):(ind(L-1)-1));
        seg2 = data1((ind(L)+1):end);
    else
        seg1 = data1((ind(i-1)+1):(ind(i)-1));
        seg2 = data1((ind(i+1)+1):(ind(i+2)-1));
    end
    
    if isempty(seg1) || isempty(seg2),
        i = i+2;
        continue;
    end
    
    lseg1 = length(seg1);
    lseg2 = length(seg2);

    dif = ind(i+1)-ind(i)+1;
    
    if dif <= npi,
        if lseg1>npint
            seg1 = seg1((lseg1-npint+1):end);
        end
        if lseg2>npint
            seg2 = seg2(1:npint);
        end
        interp = polintre (seg1, seg2, dif, ord);
        data1(ind(i):ind(i+1)) = interp;
        flag1(ind(i):ind(i+1)) = 0;
    end

    i = i+2;
end
