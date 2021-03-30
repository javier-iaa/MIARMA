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
%                              armafill.m      
%                              lincorr.m       
%                              sing.m
%                              af_simp.m       
%                              polintre.m
%                              armaint.m       
%                              pred.m          
%
% Version: 0.0.2.0
%
% Changes: 
%  - Using ft_corr optionally for a final correction of the ARMA interpolation
%
% Date: 28/03/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off all

%% Some definitions
numvers = '0.0.2.0';

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
    fprintf(2, ' #    by  J.Pascual-Granado, IAA-CSIC, Spain. 2020           #\n');
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

% default values
simpl = true;
temp = false;
always_int = true;
ft_corr = false;
nuc = 1;
mseg = 1000;
pmin = 2;
pmax = 30;
qmax = 30;
repmax = 3;
npz = 36;
npi = 4;
facmax = 6;
facmin = 4;
facint = 3;
ascii_struct = false;

if isfield( strin, 'params' )
    
    % Flag to activate the FT correction of the arma interpolation
    if isfield( strin.params, 'ft_corr' )
        ft_corr = strin.params.ft_corr;
    end
        
    % Flag to decide between af_simp.m and armafill.m
    if isfield( strin.params, 'simpl')
        simpl = strin.params.simpl;
    end

    % Default values of the parameters.
    % Min. ratio between segment length and number of parameters for
    % the model. facmin must be >= 3
    if isfield( strin.params, 'facmin')
        facmin = strin.params.facmin;
    end

    % Max. ratio between segment length and number of parameters for
    % the model. facmax = 6;
    if isfield( strin.params ,'facmax')
        facmax = strin.params.facmax;
    end

    % Lower limit in gap length for the ARMA interpolation
    % (below this limit a simpler interpolation is used)
    if isfield( strin.params, 'npi')
        npi = strin.params.npi;
    end
    
    % Lower limit in data segment length for the ARMA interpolation.
    % This must be at least d*facmin and, as min(d)=min(p+q)=pmin+0,
    % npz must be at least pmin*facmin
    if isfield( strin.params, 'npz' )
        npz = strin.params.npz;
    end
    
    % Number of iterations of the gap-filling loop. It might be interesting 
    % to use a number higher than 2 since sometimes the number of gaps
    % from one iteration to the next one is the same in spite of the gaps 
    % to fill being different.
    if isfield( strin.params, 'repmax')
        repmax = strin.params.repmax;
    end

    % Range for the search of the optimal ARMA orders [pmin,pmax] 
    if isfield( strin.params, 'pmin')
        pmin = strin.params.pmin;
    end

    if isfield( strin.params, 'pmax')
        pmax = strin.params.pmax;
    end

    % The MA order is search in the range [0,qmax]
    if isfield( strin.params, 'qmax')
        qmax = strin.params.qmax;
    end
    
    % Maximum length of the segment used to calculate ARMA order
    if isfield( strin.params, 'mseg')
        mseg = strin.params.mseg;
    end    
    
    % Number of workers to used for parallelization. Default is 1 meaning
    % no parallelization.
    if isfield( strin.params, 'nuc' )
        nuc = strin.params.nuc;
    else
        nuc = 1;
    end

    % Always interpolate or not
    if isfield( strin.params, 'always_int' )
        always_int = strin.params.always_int;
    end

    % Decides wheter to save the Akaike matrix at a temp file
    if isfield( strin.params, 'temp' )
        temp = strin.params.temp;
    end
        
    % Flag to activate (or not) ascii output with parameters and other info
    if isfield(strin.params, 'ascii_struct')
        ascii_struct = strin.params.ascii_struct;
    end

    % Full name for the file containing the Akaike matrix        
    if isfield(strin.params, 'akaname')
        akaname = strin.params.akaname;
    end
    
    % Min. ratio between interpolated datapoints and the length of the segments
    if isfield(strin.params, 'facint')
        facint = strin.params.facint;
    end
        
end

% Save parameters used in the computation in output structure for transparency
strout.params.simpl = simpl;
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

% List of parameters for armafill/af_simp
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
    %         igap = indgap(flagin);
            igap(1:2) = [];
        end
        if flagin(end)~=0
            flagin(igap(end-1):end) = 0;
    %         igap = indgap(flagin);
            igap(end-1:end) = [];
        end
    end
    
    lgaps0 = length(find(flagin~=0));
    
    % Correction of the status array for small gaps
    if npi > 1
        if verbflag
            fprintf('Step 1b - Correction of status array for small gaps\n');
        end
        [datout, flagin] = lincorr(datin, flagin, igap, npi);
    end

    igap = indgap(flagin,1);
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
    igap = indgap(flagin,1);
      
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
    
    % Optimal order (p,q) for the ARMA model of seg
%     pminstr = num2str(pmin,'%03.f');
%     pmaxstr = num2str(pmax,'%03.f');
%     qmaxstr = num2str(qmax,'%03.f');
    
    % Parallel computation was implemented a long time ago and it is 
    % deprecated now. It wasn't removed since this part of the code 
    % might be ported to other languages but, due to the restrictions in 
    % the standard version of Matlab it is not interesting to support its 
    % development here any more. So be careful when using armaord_par 
    % or any other parallelized version of the code since it is old.
    if nuc == 1
        if exist( 'akaname', 'var' )
            % fileaka is just the star id or temp
%             fileaka = akaname(1:end-4);
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

l0 = length(igap);
% l1 = l0-1;

flagout = flagin;

% Number of gaps
numgap = round(l0/2);

if numgap==1
    fprintf('\n *Starting the gap-filling iterative process*\n\n');
    fprintf('Number of gaps: %d\n', numgap);
    
    if simpl,
        [datout, flagout] = af_simp( datout, flagout, aka, igap, params, j );
    else
        [datout, flagout] = armafill( datout, flagout, aka, igap, params, j );
    end
    
    igap = indgap(flagout,j+1);
    l1 = length(igap);
    numgap = round(l1/2);
    fprintf('\nNumber of gaps remaining: %d\n', numgap);
    
else
    fprintf('\n *Starting the gap-filling iterative process*\n\n');
    fprintf('Number of gaps: %d\n\n', numgap);
    while numgap>1
        
        while (rep<repmax && l0>=2)
            if simpl,
                [datout, flagout] = af_simp( datout, flagout, aka, igap, params, j );
            else
                [datout, flagout] = armafill( datout, flagout, aka, igap, params, j );
            end
            
            igap = indgap(flagout,j+1);
            
            if isempty(igap)
%                 l1 = 0;
                numgap = 0;
                fprintf('\nNumber of gaps remaining: 0\n');
                break;
            end
            
            % Correction of the status array for small data segments
            % If you want to disable this correction just set npz to zero
            if npz > 0
                flagout = sing(flagout, npz, igap);
                igap = indgap(flagout,1);
            end
            
            % Number of gaps
            l1 = length(igap);
            numgap = round(l1/2);
            fprintf('\nNumber of gaps remaining: %d\n\n', numgap);
            
            j = j + 1;
            
            % Every iteration begins from the opposite side of the series
            datout = flipud( datout );
            flagout = fliplr( flagout );
            igap = L-igap+1;
            igap = fliplr( igap );

            % Termination condition: the number of gaps is not repeated 
            % more than twice in consecutive iterations
            if l1~=l0
                rep = 0;
                l0 = l1;
            else
                rep = rep+1;
            end
        end
        
        if mod(j,2)==0
            datout = flipud( datout );
            flagout = fliplr( flagout );
            igap = L-igap+1;
            igap = fliplr( igap );
        end

%         numgap = floor(l1/2);
        % If the number of gaps is still greater than 1 it will merge some
        % of them and repeat the main loop
        if numgap>1
            [flagout, go] = gapmerge(flagout,igap);
            if go==true
                fprintf('\n *Reinicialization with gap merging*\n\n');
                flagin(flagout==-1)=-1;
                igap = indgap(flagout,j+1);
                rep = 0;
                l0 = length(igap);
%                 l1 = l0-1;
            else
                break;
            end
        end
        j = 1;
    end
end

% if mod(j,2)==0,
%     datout = datout(end:-1:1);
% end

% Recover data segments that were taken out with sing
% datout(flagin==-1) = datin(flagin==-1);
flagout(flagin==-1 | flagin==-0.5 | flagin==0.5) = 0;

if (exist('Llin','var'))
    Llin = Llin + length(find(flagout~=0));
end

%% Final iteration using one-sided extrap (if always_int is on)
% Fill gaps left previously due to any issue in armaint that set the flag 
% <go> to False.

if (always_int && numgap > 0)
    fprintf('\n *Restarting the gap-filling iterative process with one-sided extrap*\n\n');
%     fprintf('\n *Final iteration*\n\n');
    igap = indgap(flagout);
%     [datout, flagout] = lincorr(datout, flagout, igap, inf);
    while l1>0
        if j>repmax
            fprintf(2, '\nWarning: interpolation finished before all gaps could be filled.');
            fprintf(2, '\nTry different values in the parameter structure (e.g. facint).\n\n');
            ft_corr = false;
            break;
        end
        
        [datout, flagout] = af_simp( datout, flagout, aka, igap, params, j, '1s' );
        
        igap = indgap(flagout, j+1);
        if isempty(igap)          
            fprintf('Number of gaps remaining: 0\n0');
            break;
        end
        l1 = length(igap);
        numgap = round(l1/2);
        fprintf('Number of gaps remaining: %d\n', numgap);
        
        j = j + 1;
    end
end

%% FT correction of the ARMA interpolation
if ft_corr
    datout_corr = ftcorr(datout, flagin);
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
% END
