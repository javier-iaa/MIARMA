function [flag1,go] = gapmerge(flag,ind1,varargin)
% function flag1 = gapmerge(flag,ind1,numgap) modify the flags as if the
% gaps were merged in pairs. This allow new gaps to be interpolated by
% armafill.m
% The output 'go' indicates whether any gap has been merged or not.
%
% Changes
% - Ln 84 commented
% - Minor corrections and code cleaning
%
% Version: 0.1
% Last update: 23/05/2021

% Initialization
L = length(ind1);
N = length(flag);
flag1 = flag;
go = false;
i=1;

% Limit ratio between the length of data segment and the length of gaps 
if nargin==3
    S = varargin{1};
else
    S = 4; % (default value if no input is given)
end

% first data segment around the gap merged
lseg1 = ind1(1)-1;

% second data segment around the gap merged
if L==4
    lseg2 = N - ind1(4);
else
    lseg2 = ind1(5) - ind1(4) - 1;
end

% Length of the first gap to be merged
np = ind1(4) - ind1(1);
fac1 = lseg1/np;
fac2 = lseg2/np;

% Check whether the new gap can be interpolated or not
if (fac1 > S) && (fac2 > S)
%     flag1( (1+ind1(2)):(ind1(3)-1) ) = -1;
    flag1( (1+ind1(2)):(ind1(3)-1) ) = 1;    
    go = true;
    i = i+4;
else
    i = i+2;
end

if L==4
    return;
end
    
% Now the same for the rest of the gaps
while i<L-3
    seg1 = ind1(i-1)+1:ind1(i)-1;
    lseg1 = length(seg1);
    seg2 = ind1(i+3)+1:ind1(i+4)-1;
    lseg2 = length(seg2);
    np = ind1(i+3)-ind1(i);
    fac1 = lseg1/np;
    fac2 = lseg2/np;

    if (fac1 > S) && (fac2 > S),
%         flag1( (1+ind1(i+1)):(ind1(i+2)-1) ) = -1;
        flag1( (1+ind1(i+1)):(ind1(i+2)-1) ) = 1;
        go = true;
        i = i+4;
    else
        i = i+2;
    end
end

% Last gap merging
if i==L-3
    seg1 = ind1(i-1)+1:ind1(i)-1;
    lseg1 = length(seg1);
    seg2 = ind1(i+3)+1:N;
    lseg2 = length(seg2);
    np = ind1(i+3)-ind1(i);
    fac1 = lseg1/np;
    fac2 = lseg2/np;

    if (fac1 > S) && (fac2 > S)
%         flag1( (1+ind1(i+1)):(ind1(i+2)-1) ) = -1;
        flag1( (1+ind1(i+1)):(ind1(i+2)-1) ) = 1;
        go = true;
    else
%         flag1( flag~=0 )=1;
        return;
    end
else
%     flag1( flag~=0 )=1;
    return;
end
