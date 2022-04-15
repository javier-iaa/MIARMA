function [datout, flagout, ftc] = af_simp(datin, flagin, aka, ind1, params, varargin)
% function [datout,flagout] = af_simp(datin, flagin, aka, ind1, params, varargin)
% Simplified version of armafill.m where segment length is fixed.
%
% af_simp fill gaps using ARMA models as predictors for the extrapolations
%
% Input: datin - input data array
% 
%           flagin - status array. The gaps must be correctly flagged.
% 
%           aka - Akaike coefficient matrix
% 
%           ind1 - gap indexes
% 
%           params - parameter list, is an array with these elements:
%              facmin - min. ratio between segment length and number of 
%                parameters for the model
%              facmax - max. ratio[
%              npi - inf. limit in gap length for the ARMA interpolation
%                (below this limit linear interpolation is used)
%              pmin - inf. limit por the AR order
%              fc - min. ratio between segment length and number of 
%                data points to interpolate inside the gap
% 
%           varargin{1} must be the iteration number iter
%              
%             optional inputs collected in varargin:
%              - 'lastr_aka' followed by a boolean argument
%               activate last resource solution = models with a lower
%               number of coefficients than the optimal model are used when
%               there is an insufficient amount of data.
%               - '1s' allows one-sided extrapolation when available 
%                   data does not allow forward-backward interpolation,
%                  otherwise (default), two-sided interpolation is used.
%
% Output:   datout - ARMA interpolated data series
%                flagout - residual status array
%                ftc - flag for FT correction
%
% Calls:   armaint.m
%
% Version: 0.4.2
%
% Changes from the last version: 
% - Progress numbers supressed (it is more important for armaord).
% 
% Author: Javier Pascual-Granado
%
% Date: 15/04/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flag to activate the FT correction in case armaint fails
ftc = false;

L = length(datin);
l0 = length(ind1);

% convert data into column vector
datout = reshape(datin,L,1);

%% Parameters %%
facmin = params(1);
facmax = params(2);
npi = params(3);
pmin = params(4);
fc = params(5);
iter = varargin{1};

lastr_aka_flag = find( strcmp(varargin,'lastr_aka'), 1 );
if ~isempty(lastr_aka_flag)
    lastr_aka = varargin{lastr_aka_flag+1};
else
    lastr_aka = false;
end

onesd = find( strcmp(varargin,'1s'), 1 );
if isempty( onesd )
    onesd = false;
else
    onesd = true;
end

% Data points used for the polynomial fitting.
npint = 6;

% Internal parameter S is related to the efficiency. When the length of the
% gap is > S times the length of any of the segments the algorithm is no 
%longer efficient for filling such a long gap and the gap is left unfilled.
% S = 4;

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
% i = 1;
% indlast = 0;

text_iter = sprintf('Gap filling iteration %d ----        ', iter);
fprintf(text_iter);

% Here begins the gap-filling process
while ind1f>=0,
    
%% Data segments selection
    if l1==2,   % only one gap
        if ind1(1)==1           % Left edge
            seg1 = NaN;
            subi2 = (ind1(2)+1):L;
            seg2 = datout( subi2 );
        else
            if ind1(2)==L       % Right edge
                subi1 = 1:ind1(1)-1;
                seg2 = NaN;
            else
                subi1 = 1:ind1(1)-1;
                subi2 = (ind1(2)+1):L;
                seg2 = datout( subi2 );
            end
            nf = find( flagout( subi1 ) == 1, 1, 'last');
            if ~isempty(nf),
                subi1( 1:nf ) = [];
            end
            seg1 = datout( subi1 );
        end
    else
        if ind1(1)==1          % Left edge
            seg1 = NaN;
            subi2 = (ind1(2)+1):(ind1(3)-1);
            seg2 = datout( subi2 );
        elseif ind1(2)==L,    % Right edge
            seg2 = NaN;
            subi1 = 1:ind(1)-1;
            nf = find( flagout( subi1 ) == 1, 1, 'last');
            if ~isempty(nf),
                subi1( 1:nf ) = [];
            end
%             subi1 = (indlast+1):(ind1(1)-1);
            seg1 = datout( subi1 );
        else           
            subi1 = 1:ind1(1)-1;
            nf = find( flagout( subi1 ) == 1, 1, 'last');
            if ~isempty(nf),
                subi1( 1:nf ) = [];
            end
%             subi1 = (indlast+1):(ind1(1)-1);

            if length(subi1) < 3
            % If the length of seg1 is less than 3 no ARMA model can be
            % fitted so this data segment is unusable.
            % Don't confuse this with what happens with the sing 
            % algorithm, which can be used to improve the quality of 
            % interpolations.
                if onesd
                    seg1 = NaN;
                else
                    ind1(1:2)=[];
                    l1 = length(ind1);
                    ind1f = l1-2;
                    flagout(subi1) = -1;
                    continue;
                end
                
            elseif flagout(subi1)==0.5,
            % Similarly if this segment cannot be used to perform a forward
            % extrapolation, the algorithm jumps to the next gap and for
            % the next iteration to fill this one
                if onesd
                    seg1 = NaN;
                else
                    ind1(1:2)=[];
                    l1 = length(ind1);
                    ind1f = l1-2;
                    continue;
                end
                
            else
%                 nf = find( flagout( subi1 ) == 0.5, 1, 'last' );
%                 if ~isempty(nf)
%                     subi1( 1:nf ) = [];
%                 end
                seg1 = datout( subi1 );

                subi2 = (ind1(2)+1):(ind1(3)-1);                

                if length(subi2) < 3
                % If the length of seg2 is less than 3 no ARMA model can be
                % fitted so this data segment is unusable.
                % Don't confuse this with what happens with the sing 
                % algorithm, which can be used to improve the quality of 
                % interpolations.
                    if onesd
                        seg2 = NaN;
                    else
                        ind1(1:2)=[];
                        l1 = length(ind1);
                        ind1f = l1-2;
                        flagout(subi2) = -1;
                        continue;
                    end

                elseif flagout(subi2)==-0.5
                % Similarly if this segment cannot be used to perform a forward
                % extrapolation, the algorithm jumps to the next gap and for
                % the next iteration to fill this one
                    if onesd
                        seg2 = NaN;
                    else
                        ind1(1:2)=[];
                        l1 = length(ind1);
                        ind1f = l1-2;
                        continue;
                    end
                    
                else                
%                      nf = find( flagout( subi2 ) == -0.5, 1, 'last');
%                      if ~isempty(nf),
%                          subi2( 1:nf ) = [];
%                      end
                    seg2 = datout( subi2 );
                end
            end
        end      
    end
    
    lseg1 = length(seg1);
    lseg2 = length(seg2);
    
    % number of lost datapoints in the gap
    np = ind1(2) - ind1(1) + 1;
    
    % Percentage complete
%     fprintf(repmat('\b',1,5));
%     fprintf('%3.0f %%', 100*ind1(1)/L);

 %% Perform several checks over the data
    
    % NaN conditions - no nans in seg1 and seg2
    
    nancy1 = find(isnan(seg1),1, 'last');
    nancy2 = find(isnan(seg2),1, 'first');
    nnanc1 = isempty( nancy1 );
    nnanc2 = isempty( nancy2 );
    no_nan_cond = nnanc1 && nnanc2;
   
    if no_nan_cond
        
        % Checks whether the segment length is enough to interpolate np
        % data points in the gap
           if lseg1/np<fc
                if lseg2/np<fc
                    ind1(1:2)=[];
                    l1 = length(ind1);
                    ind1f = l1 - 1;
                    continue
                else
                    if onesd
                        seg1 = NaN;
                        nnanc1 = false;
                        len = length(seg2);
                    else
                        ind1(1:2)=[];
                        l1 = length(ind1);
                        ind1f = l1 - 1;
                        continue
                    end
                end
           else
                if lseg2/np<fc
                    if onesd
                        seg2 = NaN;
                        nnanc2 = false;
                        len = length(seg1);
                    else
                        ind1(1:2)=[];
                        l1 = length(ind1);
                        ind1f = l1 - 1;
                        continue
                    end
                else
                    % Truncate segments in order to have the same length
                    difl = lseg1 - lseg2;
                    if difl > 0
                        subi1 = subi1((difl+1):end);
                        seg1 = datout( subi1 );
                    elseif difl < 0
                        subi2 = subi2(1:(lseg2+difl));
                        seg2 = datout( subi2 );
                    end
                    len = length(seg1);
                    
                    % Small gaps are linearly interpolated
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
           end
    else
        if ~nnanc1
            seg1 = seg1( (nancy1+1):end );
        else
            seg2 = seg2( 1:(nancy2-1) );
        end
        
        % Truncate segments in order to have the same length
        difl = lseg1 - lseg2;
        if difl > 0
            subi1 = subi1((difl+1):end);
            seg1 = datout( subi1 );
        elseif difl < 0
             subi2 = subi2(1:(lseg2+difl));
             seg2 = datout( subi2 );
        end
        len = length(seg1);
                    
        % Small gaps are linearly interpolated
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
                    
        % If one of the segments is nan the length of the other will be the 
        % longer necessarily.
        if lseg1 > lseg2
            len = length(seg1);
        else
            len = length(seg2);
        end
    end   
          
    % Check whether the segment length is enough to fit the arma model
    d = sum(ord);
    fac = len / d;
    if (fac <= facmin)
        
        if lastr_aka
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
                if size( ord, 1)>1
                    d = sum(ord, 2);
                    [dmin, id] = min(d);
                    p = ord(id,1);
                    q = ord(id,2);
                    ord = [p q];
                    d = dmin;
                else
                    d = sum(ord);
                end
                fac = len / d;
            end
        else
            ind1(1:2)=[];
            l1 = length(ind1);
            ind1f = l1-2;
            continue;
        end
    end

    % Too long segments are reduced by facmax for efficiency
    if (fac > facmax) && (facmax*d>fc*np),
        if nnanc1  
            newi1 = 1 + floor( (fac-facmax)*d );
            subi1 = subi1(newi1:end);
            seg1 = datout( subi1 );
        end
        if nnanc2
            newi2 = len - floor( (fac-facmax)*d );
            subi2 = subi2(1:newi2);
            seg2 = datout( subi2 );
        end
    end
            
%%  Interpolation
    
    % Interpolation algorithm. go indicates whether it was possible or not
    [interp, go] = armaint(seg1, seg2, ord, np);
    
    % Finally the interpolated segment is inserted in datout
    if go
        datout(ind1(1):ind1(2)) = interp;
        flagout(ind1(1):ind1(2)) = 0;
    else
        % If armaint could not interpolate we try with the next "optimal" order
        % If, in any case, this results insufficient we could try in the future two
        % solutions: a loop to find the order that makes it works, to restrict the 
        % orders in the MA part, since this appears to be more unstable when 
        % the q is high.
        if lastr_aka
            while ~go
                aka(cp, cq) = nan;
                minaka = min( min(aka) );
                if isnan( minaka )
                    ftc = true;
                    break
                end
                [cp, cq] = find(aka == minaka);
                q = cq - 1;
                p = cp + pmin - 1;
                ord = [p q];
                [interp, go] = armaint(seg1, seg2, ord, np);
                if go
                    datout(ind1(1):ind1(2)) = interp;
                    flagout(ind1(1):ind1(2)) = 0;
                    break
                end
            end
        else
            ftc = true;
        end
    end
        
    % Recover data points that were taken out with sing
%     reco0 = find(flagin(ind1(1):ind1(2))==-1);
%     reco = ind1(1) - 1 + reco0;
% 
%     if ~isempty(reco) && go
%         % Fix local trends that might be not modeled properly introducing jumps
%         datout = locdetrend(datout, datin, reco, reco0, seg1, seg2, interp, npint, np, ind1);
%         datout(reco) = datin(reco);
%     end
   
%     indlast = ind1(2);
    ind1(1:2)=[];
    
    % reinicialization
    l1 = length(ind1);
    ind1f = l1-2;
    aka = aka1;
    ord = ord1;
 end
 
% fprintf(repmat('\b',1,5));
fprintf('\n');

end
