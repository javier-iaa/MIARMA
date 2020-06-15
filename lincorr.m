function [data1, flag1] = lincorr(data, flag, ind, npi)
% function [data1,flag1] = lincorr(data,flag,ind,npi) use a linear
% interpolation for gap-filling the smallest gaps. 
% Inputs:   data - time series
%           flag - status array
%           ind - gap indexes
%           npi - limits the size of the gaps to be interpolated.
% Ouputs:   data1 - time series after lin.interpolation
%           flag1 - new status array
% Version: 1.0.6
% Changes from the last version: waitbars removed.
%
% Calls: polintre.m
% Author(s): Javier Pascual-Granado
% Date: 27/05/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Order used for the polynomial interpolation
ord = 3;

% Data points used for the polynomial fitting. Be careful with this, just
% a few data points is recommended.
npint = ord*2;

flag1 = flag;
data1 = data;

N = length(data);

% Change column into row
data1 = reshape(data1,1,N);

if isempty(ind)==1,
    return;
end

L = length(ind);

% h = waitbar(0,'Step 1 - Small gaps correction...');

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
%         waitbar(i/L,h);
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
%         if i==1,
%             lseg = lseg1 + dif + lseg2;            
%             interp = interp1([1:lseg1 (lseg1+dif+1):lseg],...
%                [seg1 seg2],lseg1+(1:dif),'spline');
%             interp = interp + rstd*randn(size(interp));
            
%             data1(ind(1):ind(2)) = interp;
%             flag1(ind(1):ind(2)) = 0;
%         else              
%             lseg = lseg1+dif+lseg2;            
%             interp = interp1([1:lseg1 (lseg1+dif+1):lseg],...
%                [seg1 seg2],lseg1+(1:dif),'spline');
%             interp = interp + rstd*randn(size(interp));
%             
%             data1(ind(i):ind(i+1)) = interp;
%             flag1(ind(i):ind(i+1)) = 0;
%         end
    end

    i = i+2;
%     waitbar(i/L,h);
end

% close(h);
