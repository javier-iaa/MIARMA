function [datout,flagout] = af_simp(datin, flagin, aka, ind1, params, varargin)
% function [datout,flagout] = af_simp(datin, flagin, aka, ind1, params, varargin)
% Simplified version of armafill.m where segment length is fixed.
%
% af_simp fill gaps using ARMA models as predictors for the extrapolations
%
% Input:    datin - input data array
%           flagin - status array. The gaps must be correctly flagged.
%           aka - Akaike coefficient matrix
%           ind1 - gap indexes
%           params - parameter list, is an array of 3 elements:
%              facmin - min. ratio between segment length and number of 
%                parameters for the model
%              facmax - max. ratio[
%              npi - inf. limit in gap length for the ARMA interpolation
%                (below this limit linear interpolation is used)
%              pmin - inf. limit por the AR order
%              repmax - maximum number of iterations
%
% Output:   datout - ARMA interpolated data series
%           flagout - residual status array
%
% Calls:   armaint.m
%
% Version: 0.1.9
%
% Changes from the last version: 
% - Minor fixes.
% 
% Author: Javier Pascual-Granado
%
% Date: 24/06/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = length(datin);
l0 = length(ind1);

% convert data into column vector
datout = reshape(datin,L,1);

%% Parameters %%
facmin = params(1);
facmax = params(2);
npi = params(3);
pmin = params(4);
iter = varargin{1};

% Data points used for the polynomial fitting.
npint = 6;

% Internal parameter S is related to the efficiency. When the length of the
% gap is > S times the length of any of the segments the algorithm is no 
%longer efficient for filling such a long gap and the gap is left unfilled.
S = 4;

% Akaike matrix is used to select the optimum (p,q) order. The lowest p is
% prioritized
minaka = min(min(aka));

[cp,cq] = find(aka==minaka);
q = cq - 1;
p = cp + pmin - 1;
ord = [p q];

% Static variables
aka1 = aka;
ord1 = ord;
flagout = flagin;
l1 = l0;

%% Gap-filling process %% 
ind1f = l0-1;
i = 1;

if nargin==5,
    text_iter = 'Gap filling iteration 1\n     Please wait...\n';
    fprintf(text_iter);
    fprintf('Start');
else
    text_iter = sprintf('Gap filling iteration %d\n     Please wait...\n',...
        iter);
    fprintf(text_iter);
    fprintf('Start');
end

% Here begins the gap-filling process
while ind1f>=0,
    
    % Data segments selection
    if l1==2,   % only one gap
        if ind1(i)==1,
            subi2 = (ind1(i+1)+1):L;
            seg2 = datout( subi2 );
            seg1 = NaN;
        else
            if ind1(i+1)==L
                seg2 = NaN;
            else
                subi2 = (ind1(i+1)+1):L;
                seg2 = datout( subi2 );
            end
            
            subi1 = 1:ind1(i)-1;
            nf = find( flagout( subi1 ) == 1, 1, 'last');
            if ~isempty(nf),
                subi1( 1:nf ) = [];
            end
            seg1 = datout( subi1 );
        end
    else
        if ind1(i)==1,          % Left edge
            subi2 = (ind1(i+1)+1):(ind1(i+2)-1);
            seg2 = datout( subi2 );
            seg1 = NaN;
        elseif ind1(i+1)==L,    % Right edge
            seg2 = NaN;
            
            subi1 = (ind1(i-1)+1):(ind1(i)-1);
            nf = find( flagout( subi1 ) == 1, 1, 'last');
            if ~isempty(nf),
                subi1( 1:nf ) = [];
            end
            seg1 = datout( subi1 );
        else           
            subi1 = 1:ind1(i)-1;
            % If the length of seg1 is less than 3 no ARMA model can be
            % fitted so this data segment is unusable.
            % Don't confuse this with what happens with the sing 
            % algorithm, which can be used to improve the quality of 
            % interpolations.
            if length(subi1) < 3,
                ind1(1:2)=[];
                l1 = length(ind1);
                ind1f = l1-2;
                flagout(subi1) = -1;
                continue;
            end
            % Similarly if this segment cannot be used to perform a forward
            % extrapolation, the algorithm jumps to the next gap and for
            % the next iteration to fill this one
            if flagout(subi1)==0.5,
                ind1(1:2)=[];
                l1 = length(ind1);
                ind1f = l1-2;
                continue;
            end
            nf = find( flagout( subi1 ) == 1, 1, 'last');
            if ~isempty(nf),
                subi1( 1:nf ) = [];
            end
            seg1 = datout( subi1 );

            if i==(l0-1),       % Last gap
                subi2 = (ind1(i+1)+1):L;
                seg2 = datout( subi2 );
            else
                subi2 = (ind1(i+1)+1):(ind1(i+2)-1);                
                % If the length of seg2 is less than 3 no ARMA model can be
                % fitted so this data segment is unusable.
                % Don't confuse this with what happens with the sing 
                % algorithm, which can be used to improve the quality of 
                % interpolations.

                if length(subi2) < 3,
                    ind1(1:2)=[];
                    l1 = length(ind1);
                    ind1f = l1-2;
                    flagout(subi2) = -1;
                    continue;
                end
                % Similarly if this segment cannot be used to perform a forward
                % extrapolation, the algorithm jumps to the next gap and for
                % the next iteration to fill this one
                if flagout(subi2)==-0.5,
                    ind1(1:2)=[];
                    l1 = length(ind1);
                    ind1f = l1-2;
                    continue;
                end
                
                seg2 = datout( subi2 );
            end
        end      
    end 
    
    % number of lost datapoints in the gap
    np = ind1(i+1)-ind1(i)+1;
    
    fprintf(repmat('\b',1,5));
    fprintf('%3.0f %%', 100*ind1(i)/L);


    lseg1 = length(seg1);
    lseg2 = length(seg2);
    difl = lseg1 - lseg2;

    % Data segments are truncated in order to have the same length
    if difl > 0
        subi1 = subi1((difl+1):end);
        seg1 = datout( subi1 );
    elseif difl < 0
        subi2 = subi2(1:(lseg2+difl));
        seg2 = datout( subi2 );
    end
    len = length(seg1);
    
    % Small gaps are linearly interpolated
    if isempty(find(isnan(seg1),1)) && isempty(find(isnan(seg2),1)),        
        if (np <= npi),
            if len>npint
                seg1 = seg1((len-npint+1):end);
                seg2 = seg2(1:npint);
            end
            
            interp = polintre (seg1, seg2, np, 3);
                 
            datout(ind1(1):ind1(2)) = interp;
            flagout(ind1(1):ind1(2)) = 0;
            ind1(1:2)=[];
            
            l1 = length(ind1);
            ind1f = l1-1;
            
            continue;
        end
    end
    
    % Check whether the segment length is enough to fit the arma model
    d = sum(ord);
    fac = len / d;
    if (fac <= facmin)
        [akar, akac] = find(aka);
        
        % Condition for the orders to be able to interpolate
        akaind = (akar+akac) < floor(len/facmin);
        akared = aka( akaind );
        
        if isempty(akared)
            % try the simplest ARMA model when there is only one gap
            if ind1f == 2
                p = 2; 
                q = 0;
                ord = [p q];
                d = sum(ord);
                fac = len / d;
            else
                flagout(ind1(1):ind1(2)) = 1;
                ind1(1:2)=[];
                l1 = length(ind1);
                ind1f = l1-2;

                % reinicialization
                ord = ord1;
                aka = aka1;
                continue;
            end
        else
            minaka = min(akared);
            [cp, cq] = find(aka == minaka);
            q = cq - 1;
            p = cp + pmin - 1;
            ord = [p q];
            d = sum(ord);
            fac = len / d;
        end
    end

    % Too long segments are reduced by facmax for efficiency
    if (fac > facmax) && (floor(facmax*d)>np/S),
            if isempty(find(isnan(seg1),1)),  
                newi1 = 1 + floor(fac-facmax)*d;
                subi1 = subi1(newi1:end);
                seg1 = datout( subi1 );
            end
            if isempty(find(isnan(seg2),1)),          
                newi2 = len - newi1 + 1;
                subi2 = subi2(1:newi2);
                seg2 = datout( subi2 );
            end
            len = length(seg1);
    end
                
    % Interpolation algorithm. go indicates whether it was possible or not
    [interp, go] = armaint(seg1, seg2, ord, np);
    
    % Finally the interpolated segment is inserted in datout
    if go
        datout(ind1(1):ind1(2)) = interp;
        flagout(ind1(1):ind1(2)) = 0;
    end
        
    % Recover data points that were taken out with sing
    reco0 = find(flagin(ind1(1):ind1(2))==-1);
    reco = ind1(1) - 1 + reco0;

    if ~isempty(reco) && go
        % Fix local trends that might be not modeled properly introducing jumps
        datout = locdetrend(datout, datin, reco, reco0, seg1, seg2, interp, npint, np, ind1);
        datout(reco) = datin(reco);
    end
    
    ind1(1:2)=[];
    
    % reinicialization
    l1 = length(ind1);
    ind1f = l1-2;
    aka = aka1;
    ord = ord1;
 end
 
fprintf(repmat('\b',1,5));
fprintf('\n');

end
