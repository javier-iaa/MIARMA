function flag1 = sing(flag,npz)
% function [data1,flag1] = lincorr(data,flag,ind,npi) use a linear
% interpolation for gap-filling the smallest gaps. 
% Inputs:   flag - status array
%           ind - gap indexes
%           npz - limits the size of the segments that can be used.
% Ouputs:   flag1 - new status array

% Changes from the last version:
%
%  Author(s): Javier Pascual-Granado
%  $Date: 24/02/2014 23:10$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag1 = flag;
L = length(flag);

h = waitbar(0,'Step 2- Small segments correction...');

% First data segment
indf = find(flag~=0,1);
lseg = indf-1;
if (lseg <= npz),
    flag1(1:lseg) = -1;
end
i = indf;

% First gap
indf = find(flag(i:end)==0,1);
lgap = indf-1;
i = i + lgap;

while i<L,
    waitbar(i/L,h);
    
    % Next data segment
    indf = find(flag(i:end) > 0,1);
    if isempty(indf)==0,
        lseg = indf - 1;
        if (lseg <= npz),
            flag1(i:(i+lseg-1)) = -1;
        end
        i = i + lseg;
    else
        lseg = L-i;
        if (lseg <= npz),
            flag1(i:end) = -1;
        end
        break;
    end
    
    % Next gap
    indf = find(flag(i:end)==0,1);
    lgap = indf - 1;
    i = i + lgap;
end

close(h);