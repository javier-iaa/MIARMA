function strout = MIARMA(strin)
% function strout = MIARMA(strin) 
% interpolates datapoints in a gapped time series using ARMA models to
% predict the segments of data.
% Inputs:   
%           MIARMA( strin )
%            where strin is a struct that contains the necessary inputs:
%               time, data, stat
%
%            and a set of optional inputs:
%               igap - is the gap indexes array
%               aka - is the Akaike coefficient matrix
%               temp - boolean, 1 to save temp files 0 otherwise
%
%            and parameters in the struct params that contains the fields:
%                temp, facmin, facmax, npi, npz, repmax, pmin, pmax, mseg, 
%                nuc, always_int and qmax. 
%
%            Consult the documentation for a detailed description of each 
%            of these parameters.
%
% Outputs:  If no output variable is given an ASCII file is created 
%           containing the output data and the corresponding time and
%           status
%
%           The output structure contains timeout, datout, statout
%           and also the following variables:
%
%              aka - Akaike coefficient matrix
%              ord - optimal order (p,q)
%              igap - gap indexes
%              params - a complete list of parameters used for computation
%
% Note:  All gaps must be correctly flagged for the gap-filling algorithm 
% to give an adequate output.
%
% By Javier Pascual-Granado
% <a href="matlab:web http://www.iaa.es;">IAA-CSIC, Spain</a>
%
% Dependencies:      armaord_par.m   
%                              armaord.m       
%                              indgap.m        
%                              lincorr.m       
%                              sing.m
%                              af_simp.m       
%                              polintre.m
%                              armaint.m       
%                              pred.m          
%
% Version: 0.1.0.0 
%
% Changes: 
% 
% - FT correction is activated by af_simp.m flag ftc when something 
% goes wrong.
%
% - A new flag 'lastr_aka' is introduced in order to give a higher priority
% to the optimal order (see af_simp.m for more details). This flag is 
% internal, meaning that it cannot be given as input (at least by now).
%
% - A new flag 'reco' is introduced to decide whether the datapoints 
% excluded during the gap-filling process will be recovered at the end.
%
% - Parameter simpl is removed since only af_simp.m is used now.
%
% - Merging is used for last solution resource too.
%
% - Fix: original data excluded during interpolation is now recovered.
% 
% - Multiple minor fixes and improvements, e.g. gap merging in extrap loop
%
% Date: 06/11/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Some definitions
numvers = '0.1.0.0';

lgaps0 = NaN;
Llin = NaN;

% Flag that controls screen output. Presently there are only two modes:
% 'full' and 'simple'. In the future a 'minimal' mode will be implemented 
% in order to suppress all output in screen.
if isfield( strin.params, 'verbose')
    verbose = strin.params.verbose;
else
    verbose = 'full';
end

verbflag = strcmp(verbose, 'full');

if verbflag

    % Header
    fprintf(2, '\n #############################################################\n');
    fprintf(2, ' #                                                           #\n');
    fprintf(2, ' #                    MIARMA  %s                        #\n', numvers);
    fprintf(2, ' #    by  J.Pascual-Granado, IAA-CSIC, Spain. 2021           #\n');
    fprintf(2, ' #               License GNU GPL v3.0                        #\n');
    fprintf(2, ' #                                                           #\n');
    fprintf(2, ' #############################################################\n\n');

end
    
%% Input data
timein = strin.time;
datin = strin.data;
flagin = strin.stat;

L = length(datin);

% Output variables
datout = datin;
timeout = timein;

%% Parameters

% --- Default values for parameters if no input is given ---

% Maximum length of the segment used to calculate ARMA order
mseg = 1000;

% Max. ratio between segment length and number of parameters for the model
facmax = 6;

% Min. ratio between segment length and number of parameters for the model
facmin = 4;

% Min. ratio between interpolated datapoints and the length of the segments
facint = 3;

% Number of iterations of the gap-filling loop
repmax = 3;

% Lower limit in data segment length for the ARMA interpolation.
npz = 36;

% Lower limit in gap length in order to use ARMA interpolation, below this 
% limit a simpler interpolation is used
npi = 4;

% Decides whether to save the Akaike matrix at a temp file
temp = false;

% Always interpolate or not
always_int = true;

% Flag to activate the FT correction of the arma interpolation
ft_corr = false;

% Number of workers to used for parallelization
nuc = 1;

% Range for the search of the optimal ARMA orders [pmin,pmax] 
pmin = 2;
pmax = 30;

% The MA order is search in the range [0,qmax]
qmax = 30;

% Flag to activate (or not) ascii output with parameters and other info
ascii_struct = false;

% Flag to activate recovery of excluded data points at the end of the 
% gap-filling process. False means that some data points will be substituted 
% by interpolated data points in the resulting array.
reco_flag = false;

% --- Input structure that changes parameter values ---
if isfield( strin, 'params' )
    
    if isfield( strin.params, 'ft_corr' )
        ft_corr = strin.params.ft_corr;
    end
        
    %  facmin must be >= 3
    if isfield( strin.params, 'facmin')
        facmin = strin.params.facmin;
    end

    if isfield( strin.params ,'facmax')
        facmax = strin.params.facmax;
    end

    if isfield( strin.params, 'npi')
        npi = strin.params.npi;
    end
    
    % This must be at least d*facmin and, as min(d)=min(p+q)=pmin+0,
    % npz must be at least pmin*facmin
    if isfield( strin.params, 'npz' )
        npz = strin.params.npz;
    end
    
    % It might be interesting to use a number higher than 2 since sometimes 
    % the number of gaps from one iteration to the next one is the same in 
    % spite of the gaps to fill being different.
    if isfield( strin.params, 'repmax')
        repmax = strin.params.repmax;
    end

    if isfield( strin.params, 'pmin')
        pmin = strin.params.pmin;
    end

    if isfield( strin.params, 'pmax')
        pmax = strin.params.pmax;
    end

    if isfield( strin.params, 'qmax')
        qmax = strin.params.qmax;
    end
    
    if isfield( strin.params, 'mseg')
        mseg = strin.params.mseg;
    end    
    
    if isfield( strin.params, 'nuc' )
        nuc = strin.params.nuc;
    else
        nuc = 1;
    end

    if isfield( strin.params, 'always_int' )
        always_int = strin.params.always_int;
    end

    if isfield( strin.params, 'temp' )
        temp = strin.params.temp;
    end
        
    if isfield(strin.params, 'ascii_struct')
        ascii_struct = strin.params.ascii_struct;
    end

    % Full name for the file containing the Akaike matrix      
    if isfield(strin.params, 'akaname')
        akaname = strin.params.akaname;
    end
    
    if isfield(strin.params, 'facint')
        facint = strin.params.facint;
    end
    
    if isfield(strin.params, 'reco')
        reco_flag = strin.params.reco;
    end
        
end

% Save parameters used in the computation in output structure for transparency
strout.params.temp = temp;
strout.params.always_int = always_int;
strout.params.nuc = nuc;
strout.params.mseg = mseg;
strout.params.pmin = pmin;
strout.params.pmax = pmax;
strout.params.qmax = qmax;
strout.params.repmax = repmax;
strout.params.npz = npz;
strout.params.npi= npi;
strout.params.facmax = facmax;
strout.params.facmin = facmin;
strout.params.ascii_struct = ascii_struct;
strout.params.facint = facint;
% strout.params.akaname = akaname;

% List of parameters for af_simp
params = [facmin facmax npi pmin facint];

%% Building the gap indexes
if isfield(strin, 'igap')
    % Indexes given as input
    igap = strin.igap;

else        
    % Gap indexes are calculated
    if verbflag
        fprintf('Step 1 - Finding gap indexes\n');
    end
    igap = indgap(flagin);
    
    if strcmp(always_int, false)
    % Gaps at the edges of the time series are eliminated
        if flagin(1)~=0
            flagin(1:igap(2)) = 0;
            igap(1:2) = [];
        end
        if flagin(end)~=0
            flagin(igap(end-1):end) = 0;
            igap(end-1:end) = [];
        end
    end
    
    lgaps0 = length(find(flagin~=0));
    
    % Correction of the status array for small gaps
    if npi > 1
        if verbflag
            fprintf('Step 1b - Correction for small gaps\n');
        end
        [ datout, flaglin ] = lincorr( datin, flagin, igap, npi );
    end

    flagin = flaglin;
    igap = indgap(flagin);
    if isempty(igap)
        timeout = timein;
        flagout = flagin;
        strout.timeout = timeout;
        strout.datout = datout;
        strout.statout = flagout;
        strout.igap = igap;
        return
    end
    
    % Number of linearly interpolated datapoints
    lgaps = length(find(flagin~=0));
    Llin = lgaps0 - lgaps;
    
    % Correction of the status array for small data segments
    % If you want to disable this correction just set npz to zero
    if npz > 0
        if verbflag
            fprintf('Step 2 - Small segments correction...\n');
        end
        flagin = sing(flagin, npz, igap);
    end

    % Index rebuilding
    if verbflag
        fprintf('Step 3 - Index rebuilding...\n');
    end
    igap = indgap(flagin);
      
    % If no gaps are found the program returns with no further calculations
    if isempty(igap)
        timeout = timein;
        flagout = flagin;

        if ascii_struct
            Llin = length(find(flagout~=0));
            fout = 'output.agfs';
            fich = fopen(fout,'w');
    
            fprintf(fich, '# code: MIARMA\n');
            fprintf(fich, '# version: %s\n', numvers);
            fprintf(fich, '# model (p,q): not calculated\n');
            fprintf(fich, '# length: %d\n', L);
            fprintf(fich, '# gaps_arma: 0\n');
            fprintf(fich, '# gaps_linear: %d\n', Llin);
            fprintf(fich, '# facmin: %d\n', facmin);
            fprintf(fich, '# facmax: %d\n', facmax);
            fprintf(fich, '# facint: %d\n', facint);
            fprintf(fich, '# npi: %d\n', npi);
            fprintf(fich, '# npz: %d\n', npz);
            fprintf(fich, '# niter: %d\n', repmax);
            fprintf(fich, '# mseg: %d\n', mseg);
            fprintf(fich, 'x y z\n');
   
            for i=1:L
                fprintf(fich,'%16.12f %16.13f %d\n',...
                    timeout(i), datout(i), flagout(i));
            end
            fclose(fich);
        else
            strout.timeout = timeout;
            strout.datout = datout;
            strout.statout = flagout;
            strout.igap = igap;
        end
        return;
    end
end

strout.igap = igap;

%% Search for the optimal order (p,q)
if isfield( strin, 'aka')
    % Akaike coefficient matrix given as input
    aka = strin.aka;
    
else
    % This gives the length of the largest segment without gaps
    segl = [igap(1)-1 (igap(3:2:end-1)-igap(2:2:end-2)-1) L-igap(end)];
    ML = max(segl);

    % The index of the corresponding segment
    I = find(segl==ML,1);

    % The largest segment
    if I==1
        seg = datout(1:(igap(1)-1));
    elseif I==length(segl)
        seg = datout((igap(end)+1):end);
    else
        seg = datout((igap(I*2-2)+1):(igap(I*2-1)-1));
    end
    
    % mseg is the max. length of the segment use to calculate the order
    if ML>mseg
        seg = seg(1:mseg);
    end
    
    if verbflag
        fprintf('Step 4 - Order estimation\n\nPlease wait...\n');
    end
       
    % Parallel computation was implemented a long time ago and it is 
    % deprecated now. It wasn't removed since this part of the code 
    % might be ported to other languages but, due to the restrictions in 
    % the standard version of Matlab it is not interesting to support its 
    % development here any more. So be careful when using armaord_par 
    % or any other parallelized version of the code since it is old.
    if nuc == 1
        if exist( 'akaname', 'var' )
            if verbflag
                aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
                    'qmax', qmax, 'w', akaname);
            else
                aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
                'qmax', qmax, 'verbose', 0, 'w', akaname);
            end
        
        % Note that, armaord requires the flag 'w' is the last one used
        elseif temp
            if verbflag
                aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
                    'qmax', qmax, 'w' );
            else
                aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
                'qmax', qmax, 'verbose', 0, 'w' );
            end
        else
            if verbflag
                aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
                    'qmax', qmax);
            else
                aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
                'qmax', qmax, 'verbose', 0);
            end
        end    
    else
        % Open a pool for parallel computation with nuc workers.       
        matlabpool(nuc);
        aka = armaord_par(seg, 'pmin', pmin, 'pmax', pmax);
    end
end

% The optimal (p,q) pair is found. In case of coincidence the lower p is
% the preference
[cp, cq] = find( aka == min( min( aka ) ) );
q = cq - 1;
p = cp + pmin - 1;
fprintf('\nOptimal order: [%d %d]\n\n', p, q); 

% Outputs: gap indexes and Akaike coefficient matrix
if ~ascii_struct
    strout.aka = aka;
    strout.igap = igap;
    strout.timeout = timeout;
    strout.datout = datout;
    strout.statout = flagin;
    strout.ord = [p q];
end

%% ARMA filling iterations

j = 1; % iteration-number
rep = 0; % used for the termination condition

l0 = length( igap );

flagout = flagin;
datout_tmp = datout;

% Number of gaps
numgap = l0/2;

if numgap==1
    fprintf('\n *Starting the gap-filling iterative process*\n\n');
    fprintf('Number of gaps: %d\n', numgap);
    
    [datout, flagout] = af_simp( datout, flagout, aka, igap, params,1);
    
    igap = indgap(flagout);
    l1 = length( igap );
    numgap = l1/2;
    fprintf('\nNumber of gaps remaining: %d\n', numgap);
    
else
    fprintf('\n *Starting the gap-filling iterative process*\n\n');
    fprintf('Number of gaps: %d\n\n', numgap);
    while numgap>1
        
        while (rep<repmax && l0>=2)            
            [datout, flagout, ftc] = af_simp( datout, flagout, aka, igap, ...
                params, j );
            
            % Activate the FT correction with ftc flag from af_simp
            if ftc,     ft_corr = ftc;      end
            
            igap = indgap( flagout );
            
            if isempty( igap )
                numgap = 0;
                fprintf('\nNumber of gaps remaining: 0\n');
                break;
            end
                       
            % Number of gaps
            l1 = length( igap );
            numgap = l1/2;
            fprintf('\nNumber of gaps remaining: %d\n\n', numgap);
            
            j = j + 1;

            % Termination condition: the number of gaps is not repeated 
            % more than twice in consecutive iterations
            if l1==l0
                rep = rep + 1;
            else
                rep = 0;
                l0 = l1;
            end
            
            % Every iteration begins from the opposite side of the series
            if (rep<repmax && l0>=2)
                datout = flipud( datout );
                flagout = fliplr( flagout );
                igap = L-igap+1;
                igap = fliplr( igap );
            end

        end
        
        if mod(j,2)==0
            datout = flipud( datout );
            flagout = fliplr( flagout );
            igap = L-igap+1;
            igap = fliplr( igap );
        end
        
        % If the number of gaps is still greater than 1 it will merge some
        % of them and repeat the main loop
        j = 1;
        rep = 0;
        if numgap>1
%             flagout( flagout~=1 ) = 0;
            fprintf( '\n *Reinicialization with gap merging*\n' );
            numgap0 = numgap;
            [flagout, go] = gapmerge( flagout, igap, facint );
            if go==true
%                 flagin( flagout==-1 ) = -1;
                igap = indgap( flagout );
                l0 = length( igap );
                numgap = l0/2;
                fprintf('\n Merged gaps: %d \n', numgap0-numgap );
                fprintf('\nNumber of gaps remaining: %d\n\n', numgap);
            else
                fprintf(2,'\nMerging is not effective to fill more gaps with these parameters.\n');
                break;
            end
        end
    end
end

% Recover data segments that were taken out with sing
% datout(flagin==-1) = datin(flagin==-1);
flagout( flagin~=1 ) = 0;

% if (exist('Llin','var'))
%     Llin = Llin + length(find(flagout~=0));
% end

%% Last resource solution for filling gaps 
% The optimal order condition is relaxed

if numgap==1
    fprintf( '\n *Reducing ARMA(p,q) order for the remaining gaps*\n\n' );
    [datout, flagout, ~] = af_simp( datout, flagout, aka, igap, ...
                params,1, 'lastr_aka', true );
    
    igap = indgap( flagout );
    l1 = length( igap );
    numgap = l1/2;
    fprintf('\nNumber of gaps remaining: %d\n', numgap);
    
else
    if numgap>1
        fprintf( '\n *Reducing ARMA(p,q) order for the remaining gaps*\n\n' );
    end
    
    while numgap>1
        
        while (rep<repmax && l0>=2)
            [datout, flagout, ~] = af_simp( datout, flagout, aka, igap, ...
                params, j, 'lastr_aka', true );
            
            % Activate the FT correction with ftc flag from af_simp
            if ftc,     ft_corr = ftc;      end
            
            igap = indgap( flagout );
            
            if isempty( igap )
                numgap = 0;
                fprintf('\nNumber of gaps remaining: 0\n');
                break;
            end
                       
            % Number of gaps
            l1 = length(igap);
            numgap = l1/2;
            fprintf('\nNumber of gaps remaining: %d\n\n', numgap);
            
            j = j + 1;
            
            % Every iteration begins from the opposite side of the series
            if (rep<repmax && l0>=2)
                datout = flipud( datout );
                flagout = fliplr( flagout );
                igap = L-igap+1;
                igap = fliplr( igap );
            end

            % Termination condition: the number of gaps is not repeated 
            % more than twice in consecutive iterations
            if l1~=l0
                rep = 0;
                l0 = l1;
            else
                rep = rep+1;
            end
        end
        
        % If the number of gaps is still greater than 1 it will merge some
        % of them and repeat the main loop
        j = 1; 
        rep=0;
        if numgap>1
%             flagout( flagout~=1 ) = 0;
            fprintf( '\n *Reinicialization with gap merging*\n' );
            numgap0 = numgap;
            [flagout, go] = gapmerge( flagout, igap, facint );
            if go==true
%                 flagin( flagout==-1 ) = -1;
                igap = indgap( flagout );
                l0 = length( igap );
                numgap = l0/2;
                fprintf('\n Merged gaps: %d\n', numgap0-numgap );
                fprintf('\nNumber of gaps remaining: %d\n\n', numgap);
            else
                fprintf(2,'\nMerging is not effective to fill more gaps with these parameters.\n');
                break;
            end
        end
    end
end

%% Final iteration using one-sided extrap (if always_int is on)
% Fill gaps left previously due to any issue in armaint that set the flag 
% <go> to False.

if (always_int && numgap > 0)
    fprintf('\n *Restarting the gap-filling with one-sided extrap*\n\n');
%     fprintf('\n *Final iteration*\n\n');
    igap = indgap(flagout);
%     [datout, flagout] = lincorr(datout, flagout, igap, inf);
    while l1>0
        if j>repmax
            j = 1; 
            if numgap>1
                fprintf( '\n *Reinicialization with gap merging*\n' );
                numgap0 = numgap;
                [flagout, go] = gapmerge( flagout, igap, facint );
                if go==true
                    igap = indgap( flagout );
                    l0 = length( igap );
                    numgap = l0/2;
                    fprintf('\n Merged gaps: %d\n', numgap0-numgap );
                    fprintf('\nNumber of gaps remaining: %d\n\n', numgap);
                else
                    fprintf(2,'\nMerging is not effective to fill more gaps with these parameters.\n');
                    fprintf(2, '\nWarning: interpolation finished before all gaps could be filled.');
                    fprintf(2, ['\nTry different values in the parameter structure e.g. facint, facmin, repmax ' ...
                ', also others like facmax or mseg if nonstationarity is suspected.\n\n']);
                    ft_corr = false;
                    break;
                end
            else
                fprintf(2, '\nWarning: interpolation finished before all gaps could be filled.');
                fprintf(2, ['\nTry different values in the parameter structure e.g. facint, facmin, repmax ' ...
                ', also others like facmax or mseg if nonstationarity is suspected.\n\n']);
                ft_corr = false;
                break;
            end
        end
        
        [datout, flagout] = af_simp( datout, flagout, aka, igap, params, j, '1s' );
        
        igap = indgap(flagout);
        if isempty(igap)          
            fprintf('\nNumber of gaps remaining: 0\n');
            break;
        end
        l1 = length(igap);
        numgap = l1/2;
        fprintf('\nNumber of gaps remaining: %d\n\n', numgap);
        
        j = j + 1;
    end
end

%% FT correction of the ARMA interpolation
if ft_corr
    datout_corr = ftcorr(datout, flagin);
end

% Recover original data that was excluded during the interpolation
if reco_flag
    datout(flaglin==0) = datout_tmp(flaglin==0);
    if ft_corr
        datout_corr(flaglin==0) = datout_tmp(flaglin==0);
    end
end

%% Save output
if ~ascii_struct
    % This occurs only when MIARMA is called for test purposes only
    
    if ft_corr
        strout.timeout = timeout;
        strout.datout = datout;
        strout.datout_corr = datout_corr;
        strout.statout = flagout;
    else
        strout.timeout = timeout;
        strout.datout = datout;
        strout.statout = flagout;
    end
    
else
    % This is the situation in a standard call
    fout = 'output.agfs';
    fich = fopen(fout,'w');
    
    fprintf(fich,'# code: MIARMA\n');
    fprintf(fich,'# version: %s\n', numvers);
    fprintf(fich,'# model (p,q): %d %d\n',p,q);
    fprintf(fich,'# length: %d\n',L);
    fprintf(fich,'# gaps_arma: %d\n',lgaps0-Llin);
    fprintf(fich,'# gaps_linear: %d\n',Llin);
    fprintf(fich,'# facmin: %d\n',facmin);
    fprintf(fich,'# facmax: %d\n',facmax);
    fprintf(fich, '# facint: %d\n', facint);
    fprintf(fich,'# npi: %d\n',npi);
    fprintf(fich,'# npz: %d\n',npz);
    fprintf(fich,'# niter: %d\n',repmax);
    fprintf(fich,'# mseg: %d\n\n',mseg);
    fprintf(fich,'x y z\n');
    
    if ft_corr
        for i=1:length(timeout)
            fprintf(fich,'%16.12f %16.13f %d\n',...
                timeout(i),datout_corr(i),flagout(i));
        end
    else       
        for i=1:length(timeout)
            fprintf(fich,'%16.12f %16.13f %d\n',...
                timeout(i),datout(i),flagout(i));
        end
    end
    fclose(fich);
    
    fprintf('\n  Interpolation finished successfully.  \n');
end

end
% END
