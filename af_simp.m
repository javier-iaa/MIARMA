function [datout, flagout, ftc] = af_simp(strin, iter, varargin)
% function [datout,flagout] = af_simp(strin, iter, varargin)
% Fill the gaps using ARMA models as predictors for the extrapolations.
%
%
% Input:    - strin - input structure as defined in MIARMA.m
%           - iter - (for debug) is a pair [i j] where i is the file number 
%            deb00i.csv and j is the folder 00j
%              
% Optional inputs collected in varargin:
%           - 'lastr_aka' followed by a boolean argument activate last 
%              resource solution = models with a lower number of 
%              coefficients than the optimal model are used when there is 
%              an insufficient amount of data.
%           - '1s' allows one-sided extrapolation when available 
%              data does not allow forward-backward interpolation,
%              otherwise (default), two-sided interpolation is used.
%           - 'nw' followed by the number of workers for the parallel loop.
%              Default: min(numcores,4) to keep the memory usage low, since
%              every worker is a full MATLAB process (0 = serial execution).
%
% Output:   - datout - ARMA interpolated data series
%           - flagout - residual status array
%           - ftc - flag for FT correction
%
% Note: af_simp.m is a variation of armafill.m where the length of the
% segments was asymmetric, now it is fixed to avoid introducing biases.
%
% Calls:   armaint.m
%          info2table.m
%
% Version: 0.5.0
%
% Changes from the last version: 
% - The gap-filling loop has been parallelized with parfor: the sequential
% <while ind1f>=0> loop was replaced by a <parfor k=1:ng> over the gaps.
% Every gap is filled independently using the original (unfilled) data, so
% the interpolated values of the already processed gaps are no longer used
% to build the segments of the following ones; results may slightly differ
% from the sequential algorithm.
% - The debug information is collected inside the parallel loop and written
% to the tables (info2table) once the loop has ended.
% - Added a guard that skips any <igap> entry whose range is not a real
% gap (flagin(g1:g2)==1). This prevents an already filled segment from
% being overwritten and re-flagged as a gap when af_simp is called again
% with a stale <igap> after the filling has converged.
% 
% Author: Javier Pascual-Granado
%
% Date: 26/09/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Flag to activate the FT correction in case armaint fails
%ftc = false;

% Unpack the structure
datin = strin.datout;
flagin = strin.statout;
ind1 = strin.igap;
facmin = strin.params.facmin;
facmax = strin.params.facmax;
npi= strin.params.npi;
pmin = strin.params.pmin;
fc = strin.params.facint;
mem = strin.params.mem;
aka = strin.aka;

L = length(datin);
l0 = length(ind1);

% convert data into column vector
datout = reshape(datin,L,1);

%% Optional inputs %%

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

% Number of workers for the parfor loop. Passed as 'nw', a nonnegative
% integer (0 = serial on the client, N = pool with N workers). If the
% existing pool has fewer workers than requested, its size is used.
nwreq = find( strcmp(varargin,'nw'), 1 );
if ~isempty(nwreq)
    nwval = varargin{nwreq+1};
else
    % Memory conscious default: never open more than 4 worker processes
    % (each one is a full MATLAB process with copies of the data).
    nwval = min(feature('numcores'), 4);
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

%% Gap-filling process %% 
text_iter = sprintf('Gap filling iteration %d.%d ----        ', iter(2), iter(1));
fprintf(text_iter);

% Number of gaps to fill
ng = l0/2;

% Boundaries of every gap: gs = start index, ge = end index and ns = start
% index of the next gap (NaN for the last one)
gs = ind1(1:2:end);
ge = ind1(2:2:end);
ns = [ ind1(3:2:end) nan ];

 % Here begins the gap-filling process (parallel version)
 %
 % NOTE: in the parfor version every gap is filled independently using the
 % original (unfilled) data as input. Contrary to the sequential version,
 % the interpolated values of the already processed gaps are not used to
 % build the segments of the following ones. This makes the computation
 % parallelizable at the cost of (slightly) shorter segments, since the NaNs
 % of the remaining gaps are truncated as usual. The result for the first
 % gap coincides with the sequential one.

% Containers collecting the results of every gap:
% - Dfills{k}: interpolated values for gap k (empty if not interpolated)
% - Fwrites{k}: flag update [r1 r2 value] applied to flagout after the loop
% - Dinfos{k}: debug info from armaint to be written after the loop
% - Ftc_gap(k): true when the FT correction is required for gap k
Dfills = cell(ng,1);
Fwrites = cell(ng,1);
Dinfos = cell(ng,1);
Ftc_gap = false(ng,1);

% In debug mode the loop is executed serially on the client (nw=0) because
% armaint may start its own pool (armax_par), which is not allowed inside
% the parfor workers. Otherwise the default parallel pool is used.
if strin.flags.debug_flag
    % In debug mode run serially on the client because armaint may start
    % its own pool (armax_par), not allowed inside parfor workers.
    nw = 0;
elseif nwval == 0
    % Explicit serial request
    nw = 0;
else
    % If a pool already exists, respect its size (the user controls it with
    % parpool). Otherwise create a pool with exactly the requested number
    % of workers: passing a number to parfor alone does not limit the pool
    % size, which would otherwise be spawned at the cluster default and
    % waste memory.
    p = gcp('nocreate');
    if isempty(p)
        nw = min(nwval, feature('numcores'));
        parpool(nw);
    else
        nw = p.NumWorkers;
    end
end

parfor (k = 1:ng, nw)

    % Number of gaps remaining (including the current one) and boundaries
    rem = ng - k + 1;
    g1 = gs(k);
    g2 = ge(k);
    if k<ng
        nn = ns(k);
    else
        nn = nan;
    end

    % Guard against stale gaps: a real gap has flagin==1 over its whole
    % range. If the range is already filled (flag 0) or otherwise not a gap
    % (e.g. igap was recomputed on a previous call and then reused), the
    % gap is simply skipped, otherwise an already filled segment could be
    % overwritten and re-flagged as a gap.
    if ~all( flagin(g1:g2)==1 )
        continue;
    end

    %seg1 = [];
    seg2 = [];
    subi1 = [];
    subi2 = [];

    % Local copies of the (re)initialized model
    ord = ord1;
    aka = aka1;
    Dinfos{k} = cell(0,2);

    % Local temporaries used in the (optional) order rescan
    cp = [];
    cq = [];

 %% Data segments selection
    if rem==1   % only one gap
        if g1==1           % Left edge
            seg1 = nan;
            subi2 = (g2+1):L;
            seg2 = datout( subi2 );
        else
            if g2==L       % Right edge
                subi1 = 1:g1-1;
                seg2 = nan;
            else
                subi1 = 1:g1-1;
                subi2 = (g2+1):L;
                seg2 = datout( subi2 );
            end
            nf = find( flagin( subi1 ) == 1, 1, 'last');
            if ~isempty(nf)
                subi1( 1:nf ) = [];
            end
            seg1 = datout( subi1 );
        end
    else
        if g1==1        % Left edge
            seg1 = nan;
            subi2 = (g2+1):(nn-1);
            seg2 = datout( subi2 );
        elseif g2==L    % Right edge
            seg2 = nan;
            subi1 = 1:g1-1;
            nf = find( flagin( subi1 ) == 1, 1, 'last');
            if ~isempty(nf)
                subi1( 1:nf ) = [];
            end
            seg1 = datout( subi1 );
        else
            subi1 = 1:g1-1;
            nf = find( flagin( subi1 ) == 1, 1, 'last');
            if ~isempty(nf)
                subi1( 1:nf ) = [];
            end

            if length(subi1)<3
            % If the length of seg1 is less than 3 no ARMA model can be
            % fitted so this data segment is unusable.
            % Don't confuse this with what happens with the sing
            % algorithm, which can be used to improve the quality of
            % interpolations.
                if onesd
                    seg1 = nan;
                else
                    if ~isempty(subi1)
                        Fwrites{k} = [subi1(1), subi1(end), -1];
                    end
                    continue;
                end

            elseif flagin(subi1)==0.5
            % Similarly if this segment cannot be used to perform a forward
            % extrapolation, the algorithm jumps to the next gap (in the
            % parallel version the gap is simply left unfilled)
                if onesd
                    seg1 = nan;
                else
                    continue;
                end

            else

                seg1 = datout( subi1 );
                subi2 = (g2+1):(nn-1);

                if length(subi2) < 3
                % If the length of seg2 is less than 3 no ARMA model can be
                % fitted so this data segment is unusable.
                % Don't confuse this with what happens with the sing
                % algorithm, which can be used to improve the quality of
                % interpolations.
                    if onesd
                        seg2 = nan;
                    else
                        if ~isempty(subi2)
                            Fwrites{k} = [subi2(1), subi2(end), -1];
                        end
                        continue;
                    end

                elseif flagin(subi2)==-0.5
                % Similarly if this segment cannot be used to perform a
                % forward extrapolation, the algorithm jumps to the next
                % gap (in the parallel version the gap is left unfilled)
                    if onesd
                        seg2 = nan;
                    else
                        continue;
                    end

                else
                    seg2 = datout( subi2 );
                end
            end
        end
    end

    lseg1 = length(seg1);
    lseg2 = length(seg2);

    % number of lost datapoints in the gap
    np = g2 - g1 + 1;

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
                    continue
                else
                    if onesd
                        seg1 = nan;
                        nnanc1 = false;
                        len = length(seg2);
                    else
                        continue
                    end
                end
           else
                if lseg2/np<fc
                    if onesd
                        seg2 = nan;
                        nnanc2 = false;
                        len = length(seg1);
                    else
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
                    if (np <= npi)
                        if len>npint
                            seg1 = seg1((len-npint+1):end);
                            seg2 = seg2(1:npint);
                        end

                        interp = polintre (seg1, seg2, np, 3);

                        Dfills{k} = interp;
                        Fwrites{k} = [g1, g2, 0];

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
        if (np <= npi)
            if len>npint
                seg1 = seg1((len-npint+1):end);
                seg2 = seg2(1:npint);
            end

            interp = polintre (seg1, seg2, np, 3);

            Dfills{k} = interp;
            Fwrites{k} = [g1, g2, 0];

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
                % try the simplest ARMA model when there are two remaining
                % gaps (the sequential code checked this through the ind1f
                % counter, which cannot be kept in the parallel version)
                if rem == 2
                    p = 2;
                    q = 0;
                    ord = [p q];
                    d = sum(ord);
                    fac = len / d;
                else
                    Fwrites{k} = [g1, g2, 1];
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
            continue;
        end
    end

    % Too long segments are reduced by facmax for efficiency
    if (fac > facmax) && (facmax*d>fc*np)
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
    % Set the debug flag for testing purposes
    if strin.flags.debug_flag
        [interp, go, info] = armaint(seg1, seg2, ord, np, 'mem', mem, ...
            'debug');
        Dinfos{k} = [ Dinfos{k}; {info, true} ];
    else
        [interp, go] = armaint(seg1, seg2, ord, np, 'mem', mem);
    end

    % Finally the interpolated segment is inserted in datout
    if go
        Dfills{k} = interp;
        Fwrites{k} = [g1, g2, 0];
    else
        % If armaint could not interpolate we try with next "optimal" order
        % If, in any case, this results insufficient we could try in the
        % future two solutions: a loop to find the order that makes it
        % works, to restrict the orders in the MA part, since this appears
        % to be more unstable when the q is high.
        if lastr_aka
            while ~go
                aka(cp, cq) = nan;
                minaka = min( min(aka) );
                if isnan( minaka )
                    Ftc_gap(k) = true;
                    break
                end
                [cp, cq] = find(aka == minaka);
                q = cq - 1;
                p = cp + pmin - 1;
                ord = [p q];
                if strin.flags.debug_flag
                    [interp, go, info] = armaint(seg1, seg2, ord, np, ...
                        'mem', mem, 'debug');
                    if go
                        Dinfos{k} = [ Dinfos{k}; {info, false} ];
                        Dfills{k} = interp;
                        Fwrites{k} = [g1, g2, 0];
                        break
                    end
                else
                    [interp, go] = armaint(seg1, seg2, ord, np, ...
                        'mem', mem);
                    if go
                        Dfills{k} = interp;
                        Fwrites{k} = [g1, g2, 0];
                        break
                    end
                end
            end
        else
            Ftc_gap(k) = true;
        end
    end
end

 %% Update datout and flagout with the results of every gap. The flag
 %  updatings are applied in the same order as in the sequential version so
 %  that overwritings produce the same final state.
for k = 1:ng
    if ~isempty(Dfills{k})
        datout( gs(k):ge(k) ) = Dfills{k};
    end
    if ~isempty(Fwrites{k})
        flagout( Fwrites{k}(1):Fwrites{k}(2) ) = Fwrites{k}(3);
    end
end

% Write the debug information sequentially to avoid concurrent writes from
% the parallel workers (info2table appends rows to csv files)
if strin.flags.debug_flag
    gnit = 1;
    for k = 1:ng
        for j = 1:size(Dinfos{k},1)
            if Dinfos{k}{j,2}
                info2table(Dinfos{k}{j,1}, gnit, iter);
            else
                info2table(Dinfos{k}{j,1}, gnit);
            end
            gnit = gnit + 1;
        end
    end
end

% FT correction required in at least one of the gaps
ftc = any( Ftc_gap );
 
% fprintf(repmat('\b',1,5));
fprintf('\n');

end
