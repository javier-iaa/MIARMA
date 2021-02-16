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
% Version: 0.0.0.2
%
% Changes: 
% - Header in red color.
% - Control the ratio between num of points fed to models and 
%  num of points interpolated through parameter facint.
% - New parameter qmax.
%
% Date: 06/02/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off all

%% Some definitions
numvers = '0.0.0.2';

lgaps0 = NaN;
Llin = NaN;

% Header
fprintf(2, '\n #############################################################\n');
fprintf(2, ' #                                                           #\n');
fprintf(2, ' #                    MIARMA  %s                        #\n', numvers);
fprintf(2, ' #    by  J.Pascual-Granado, IAA-CSIC, Spain. 2020           #\n');
fprintf(2, ' #               License GNU GPL v3.0                        #\n');
fprintf(2, ' #                                                           #\n');
fprintf(2, ' #############################################################\n\n');

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
    
    % Number of iterations of the gap-filling loop
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
    fprintf('Step 1 - Finding gap indexes\n');
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
        fprintf('Step 1b - Correction of status array for small gaps\n');
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
        fprintf('Step 2 - Small segments correction...\n');
        flagin = sing(flagin, npz, igap);
    end

    % Index rebuilding
    fprintf('Step 3 - Index rebuilding...\n');
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
    
    fprintf('Step 4 - Order estimation\n\nPlease wait...\n');
    
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
            fileaka = akaname(1:end-4);
            aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
                'qmax', qmax, 'w', fileaka);
        elseif temp            
            aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
                'qmax', qmax, 'w' );
        else            
            aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
                'qmax', qmax);            
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
fprintf('Number of gaps: %d\n', numgap);

if numgap==1
    if simpl,
        [datout, flagout] = af_simp( datout, flagout, aka, igap, params, j );
    else
        [datout, flagout] = armafill( datout, flagout, aka, igap, params, j );
    end
    
%     if mod(j,2)==0
%         datout = datout(end:-1:1);
%     end

    l1 = length(igap);
    numgap = round(l1/2);
    fprintf('Number of gaps: %d\n', numgap);
    
else
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
                fprintf('Number of gaps: 0\n');
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
            fprintf('Number of gaps: %d\n', numgap);
            
            j = j + 1;
            
            % Every iteration begins from the opposite side of the series
            datout = datout(end:-1:1);
            flagout = flagout(end:-1:1);
            igap = L-igap+1;
            igap = igap(end:-1:1);

            % Termination condition: the number of gaps is not repeated 
            % more than twice in consecutive iterations
            if l1~=l0
                rep = 0;
            else
                rep = rep+1;
            end
            l0 = l1;
        end
        
        if mod(j,2)==0
            datout = datout(end:-1:1);
            flagout = flagout(end:-1:1);
            igap = L-igap+1;
            igap = igap(end:-1:1);
        end

%         numgap = floor(l1/2);
        % If the number of gaps is still greater than 1 it will merge some
        % of them and repeat the main loop
        if numgap>1
            [flagout, go] = gapmerge(flagout,igap);
            if go==true
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

%% Here goes a new section for always_int which fills the gaps left due to 
% unsolved issues in armaint that activate the flag <go> to False
% When ARMA cannot be used we will use other algorithm (polintre)
% If the flag always_int is false this section should not be run.

igap = indgap(flagout);
[datout, flagout] = lincorr(datout, flagout, igap, inf);

%% Save output
if ~ascii_struct
    % This occurs only when MIARMA is called for test purposes only
    strout.timeout = timeout;
    strout.datout = datout;
    strout.statout = flagout;
    
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
    
    for i=1:length(timeout)
        fprintf(fich,'%16.12f %16.13f %d\n',...
            timeout(i),datout(i),flagout(i));
    end
    fclose(fich);
    
    fprintf('\n  Interpolation finished successfully.  \n');
end
% END
