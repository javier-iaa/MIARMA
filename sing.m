function flag1 = sing(flag, npz, igap)
% function flag1 = sing(flag, npz, igap) correct data segments that are too
% small to fit models that can be used to interpolate the gaps reliably.
%
% ARMA models are fitted if the data segments have at least
% gaplenfac*length_of_gap and more than npz
%
% Inputs:   flag - status array
%           npz - limits the size of the segments that can be used.
%           igap - gap indexes
% Outputs:   flag1 - new status array

% Changes from the last version:
%   - Rebuilt from scratch
%   - New criterion to fix data segments is introduced: gap length. 
%   - Several minor fixes.
%
%  Author(s): Javier Pascual-Granado
%  $Date: 25/10/2019$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% A factor of 1.0 ensures a duty cycle of 66% for a single interpolation
gaplenfac = 1.0;

flag1 = flag;
L = length(flag);
lg = length(igap);

h = waitbar(0,'Step 2- Small segments correction...');

% First data segment and gap
lseg = igap(1) - 1;
lgap = igap(2) - igap(1) + 1;
if (lseg <= npz || lseg <= gaplenfac*lgap),
    flag1(1:lseg) = 0.5;   
    % 0.5 means that the gap on the right cannot be filled
end

i = igap(2) + 1; % start of the segment

for k = 3:2:lg,
    waitbar(i/L,h);
    
    % Next data segment
    lseg = igap(k) - i;
    if (lseg <= npz || lseg <= gaplenfac*lgap),
        flag1(i:(igap(k)-1)) = -0.5;
        % 0.5 means that the gap on the left cannot be filled
    end
    lgap = igap(k+1) - igap(k) + 1;
    if (lseg <= npz || lseg <= gaplenfac*lgap),
        if flag1(i:(igap(k)-1)) == -0.5
           flag1(i:(igap(k)-1)) = 1; % if data cannot be used to fill any gap
        else
            flag1(i:(igap(k)-1)) = 0.5;
        end
    end
    
    i = igap(k+1)+1; % start of the next segment
    
end

lseg = L - i + 1;
if (lseg <= npz || lseg <= gaplenfac*lgap),
    flag1(i:end) = -0.5;
end

close(h);