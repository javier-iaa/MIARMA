function outStruct = MIARMA(inStruct)
% function outStruct = MIARMA(inStruct) 
% Interpolates datapoints in a gapped time series using ARMA models to
% predict the segments of data that are imputed.
%
% Inputs:   
%           MIARMA( inStruct )
%            where inStruct is a struct that must contains the inputs:
%               time, data, stat (not mandatory)
%
%            Optional inputs are:
%               igap - is the gap indexes array
%               aka - is the Akaike coefficient matrix
%               temp - boolean, 1 to save temp files 0 otherwise
%
%            Parameters that can be passed through inStruct.params:
%                facmin, facmax, npi, npz, pmin, pmax, qmax, mseg,
%                cutoff_level
%   
%            Flags that can be passed through inStruct.flags:
%                temp, always_int, ft_corr, verbose, ascii_struct, 
%                reco, debug
%
%           Alternatively, inStruct can be a string containing a 
%           filename of a two or three columns data file with floating 
%           number notation assumed.
%
%           An .ini file containing the parameters can be optionally passed
%           to the program. The only requirements for the ini file are
%           these:
%           - each parameter on a different line following this notation:
%           - name value
%           - all parameters should have a value (1 or 0 if logical)
%           - last line must be a carriage return
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
% Dependencies:                armaord.m       
%                              indgap.m        
%                              lincorr.m       
%                              sing.m
%                              af_simp.m       
%                              polintre.m
%                              armaint.m       
%                              pred.m
%                              autoarmaord.m
%                              fastCGSA.m
%                              saveout.m
%                              defpars.m
%
% Version: 0.1.2.10
%
% Changes:
% - BUGFIX: gap merging was activated during the 2nd ARMA filling run 
% before expected.
% - BUGFIX: igap incorrectly passed to af_simp during 3rd ARMA filling run.
% - Other minor fixes and format improvements.
%
% Date: 30/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numvers = '0.1.2.10';

%% Warning messages

% Disable all standard warnings. Specific warning messages are created below.
warning('off','all')

warning_m1 = [ '\nWarning: interpolation finished before all gaps could be filled.' ...
                '\nTry different values in the parameter structure e.g. facint, facmin, npi' ...
                ', also others like facmax or mseg if nonstationarity is suspected.\n\n'];
            
warning_m2 = '\nWarning: computing time could be up to several hours.\n\n';

warning_m3 = '\nMerging is not effective to fill more gaps with these parameters.\n';

%% Input data
if ischar( inStruct )
    filename = inStruct;
    outStruct.filename = filename;
       
    % Here data is imported from an ASCII file having 3 columns: time, flux
    % and status
    data = importdata(filename);

    % Name without extension
    akaname = filename(1:end-4);
    resFolder = akaname;

    % Depending on the characteristics of the file, importdata may
    % generate a scalar structure that is here converted into a matrix
    if isstruct(data)
        outStruct.time = data.data(:,1);
        outStruct.data = data.data(:,2);
    else
        outStruct.time = data(:,1);
        outStruct.data = data(:,2);
    end

    % Look for ini file to import parameters
    inifile = sprintf("%s.ini", akaname);
    if isfile(inifile)

        % List of available parameters
        parlist = {'mem', 'folder', 'ft_corr', 'facmin', 'facmax', 'npi', ...
            'npz', 'pmin', 'pmax', 'qmax', 'mseg', 'always_int', 'temp', ...
            'ascii_struct', 'akaname', 'facint', 'reco', 'cutoff', ...
            'debug', 'verbose'};
        
        ini = fopen(inifile, "r");
        iniln = fgetln(ini);

        while ~strcmp(iniln, "\n")
            inipar = split(iniln);
            parname = inipar{1};
            parval = inipar{2};
            if any( strcmp(parname, parlist) )
                if ~any( strcmp(parname, {'folder', 'akaname'}) )
                    parval = str2double(parval);
                end
                outStruct.params.(parname) = parval;
            end
            iniln = fgetl(ini);
        end
    end
end

tc = inStruct.time;
sc = inStruct.data;

if isfield(inStruct, 'stat')
    statc = inStruct.stat;
else
    statc = zeros( size(tc) );
end

% Sampling regularization
[timein, datin, flagin] = regsamp(tc, sc, statc);

frac_seg = maxseg(datin, flagin);

% Use fractal fraction to warn about long computing time
[~, ~, ~, fracfrac, ~] = fastCGSA( frac_seg );
if fracfrac<80
    fprintf(2, warning_m2);
    rep_lim = 10; % used for deterministic signals
else
    rep_lim = 2;
end

L = length(timein);
datout = datin;

% Populate outStruct with default parameters and flags
outStruct.params = defpars('params');
outStruct.flags = defpars('flags');

% Flag that controls screen output. Presently there are only two modes:
% 'full' and 'none'. In the future a 'minimal' mode will be implemented 
% in order to suppress most output in screen.
if isfield(inStruct, 'flags')
    if isfield(inStruct.flags, 'verbose')
        outStruct.flags.verbose = inStruct.flags.verbose;
    else
        outStruct.flags.verbose = 'none';
    end
else
    outStruct.flags.verbose = 'none';
end

verbflag = strcmp(outStruct.flags.verbose, 'full');

if verbflag

    % Header
    fprintf(2, '\n #################################################\n');
    fprintf(2, ' #                                               #\n');
    fprintf(2, ' #                 MIARMA  %s               #\n', numvers);
    fprintf(2, ' #  by J.Pascual-Granado, IAA-CSIC, Spain. 2026  #\n');
    fprintf(2, ' #              License GNU GPL v3.0             #\n');
    fprintf(2, ' #                                               #\n');
    fprintf(2, ' #################################################\n');

end
    
%% Parameters

% --- Default values for parameters if no input is given ---

% Output folder (numbered)
% Important: if there exists already a result folder with another name 
% there must be param.folder input parameter inStructg setting it
if ~ischar( inStruct )
    resList = dir('res*');
    if ~isempty(resList)
        lastFolder = resList(end).name;
        numFolder = str2double( lastFolder(4:end) );
        resFolder = sprintf('res%0.3d', numFolder+1);
    else
        resFolder = 'res001';
    end
end

% --- Input structure that changes parameter values ---
if isfield( inStruct, 'params' )

    % Total system memory (to avoid overflow issues)
    if isfield( inStruct.params, 'mem')
        outStruct.params.mem = inStruct.params.mem;
    end

    %  Set the output folder
    if isfield( inStruct.params, 'folder' )
        resFolder = inStruct.params.folder;
        outStruct.params.resFolder = resFolder;
    end

    %  facmin must be >= 3
    if isfield( inStruct.params, 'facmin')
        outStruct.params.facmin = inStruct.params.facmin;
        if inStruct.params.facmin < 3
            fprintf(2,' Warning: facmin < 3   This cannot go well!\n\n');
        end
    end

    if isfield( inStruct.params ,'facmax')
        outStruct.params.facmax = inStruct.params.facmax;
    end

    if isfield( inStruct.params, 'npi')
        outStruct.params.npi = inStruct.params.npi;
    end
    
    % This must be at least d*facmin and, as min(d)=min(p+q)=pmin+0,
    % npz must be at least pmin*facmin
    if isfield( inStruct.params, 'npz' )
        outStruct.params.npz = inStruct.params.npz;
    end
    
    if isfield( inStruct.params, 'pmin')
        pmin = inStruct.params.pmin;
        outStruct.params.pmin = pmin;
        outStruct.flags.auto_flag = false;
    else
        pmin = outStruct.params.pmin;
    end

    if isfield( inStruct.params, 'pmax')
        pmax = inStruct.params.pmax;
        outStruct.params.pmax = pmax;
        outStruct.flags.auto_flag = false;
    else
        pmax = outStruct.params.pmax;
    end

    if isfield( inStruct.params, 'qmax')
        qmax = inStruct.params.qmax;
        outStruct.params.qmax = qmax;
        outStruct.flags.auto_flag = false;
    else
        qmax = outStruct.params.qmax;
    end
    
    if isfield( inStruct.params, 'mseg')
        mseg = inStruct.params.mseg;
        outStruct.params.mseg = mseg;
    else
        mseg = outStruct.params.mseg;
    end    

    % Full name for the file containing the Akaike matrix      
    if isfield(inStruct.params, 'akaname')
        akaname = inStruct.params.akaname;
        outStruct.params.akaname = akaname;
    else
        akaname = outStruct.params.akaname;
    end
    
    if isfield(inStruct.params, 'facint')
        outStruct.params.facint = inStruct.params.facint;
    end
        
    if isfield(inStruct.params, 'cutoff')
        outStruct.params.cutoff_level = inStruct.params.cutoff;
    end
    
end

% Flags structure
if isfield( inStruct, 'flags')

    %  Fourier correction
    if isfield( inStruct.flags, 'ft_corr' )
        outStruct.flags.ft_corr = inStruct.flags.ft_corr;
    end

    if isfield( inStruct.flags, 'always_int' )
        outStruct.flags.always_int = inStruct.flags.always_int;
    end

    if isfield( inStruct.flags, 'temp' )
        outStruct.flags.temp = inStruct.flags.temp;
    end
        
    if isfield(inStruct.flags, 'ascii_struct')
        outStruct.flags.ascii_struct = inStruct.flags.ascii_struct;
    elseif exist('filename','var')
        outStruct.flags.ascii_struct = true;
    end

    if isfield(inStruct.flags, 'reco')
        outStruct.flags.reco_flag = inStruct.flags.reco;
    end

    if isfield(inStruct.flags, 'debug')
        outStruct.flags.debug_flag = inStruct.flags.debug;
    end
end   

if ~ischar( inStruct )
    if verbflag
        fprintf(2,'\n Using folder %s for results.\n\n', resFolder);
    end
end

% Save version number
outStruct.numvers = numvers;

%% Prepare folder
if ~isfield(inStruct, 'aka')
    if ~isfolder(resFolder)
        mkdir(resFolder);
    end
    cd(resFolder);
end

%% Building the gap indexes
if isfield(inStruct, 'igap')
    % Indexes given as input
    igap = inStruct.igap;

else        
    % Gap indexes are calculated (first and last inside the gap)
    if verbflag
        fprintf('Step 1 - Finding gap indexes\n');
    end
    igap = indgap(flagin);
    
    if strcmp(outStruct.flags.always_int, false)
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
    if outStruct.params.npi > 1
        if verbflag
            fprintf('Step 1b - Correction for small gaps\n');
        end
        [ datout, flaglin ] = lincorr( datin, flagin, igap, outStruct.params.npi );
    end

    flagin = flaglin;
    flagout = flagin;
    igap = indgap(flagin);

    % Number of linearly interpolated datapoints
    lgaps = length(find(flagin~=0));
    Llin = lgaps0 - lgaps;

    % Save parameters in output structure
    outStruct.lgaps0 = lgaps0;
    outStruct.Llin = Llin;
    outStruct.L = L;
    outStruct.timeout = timein;
    outStruct.datout = datout;
    outStruct.statout = flagout;
    outStruct.igap = igap;

    % End the program if no gaps are found
    if isempty(igap)
        saveout(outStruct);
        return
    end
    
    % Correction of the status array for small data segments
    % If you want to disable this correction just set npz to zero
    if outStruct.params.npz > 0
        if verbflag
            fprintf('Step 2 - Small segments correction...\n');
        end
        flagin = sing(flagin, outStruct.params.npz, igap);
        flagout = flagin;
        outStruct.statout = flagout;
    end

    % Index rebuilding
    if verbflag
        fprintf('Step 3 - Index rebuilding...\n');
    end
    igap = indgap(flagin);
    outStruct.igap = igap;
      
    % If no gaps are found the program returns with no further calculations
    if isempty(igap)
        % flagout = flagin;

        if outStruct.flags.ascii_struct
            Llin = length(find(flagout~=0));
            fout = 'output.agfs';
            fich = fopen(fout,'w');
    
            fprintf(fich, '# code: MIARMA\n');
            fprintf(fich, '# version: %s\n', numvers);
            fprintf(fich, '# model (p,q): not calculated\n');
            fprintf(fich, '# length: %d\n', L);
            fprintf(fich, '# gaps_arma: 0\n');
            fprintf(fich, '# gaps_linear: %d\n', Llin);
            fprintf(fich, '# facmin: %d\n', outStruct.params.facmin);
            fprintf(fich, '# facmax: %d\n', outStruct.params.facmax);
            fprintf(fich, '# facint: %d\n', outStruct.params.facint);
            fprintf(fich, '# npi: %d\n', outStruct.params.npi);
            fprintf(fich, '# npz: %d\n', outStruct.params.npz);
            fprintf(fich, '# mseg: %d\n', outStruct.params.mseg);
            fprintf(fich, 'x y z\n');
   
            for i=1:L
                fprintf(fich,'%16.12f %16.13f %d\n',...
                    timein(i), datout(i), flagout(i));
            end
            fclose(fich);
        % else
        %     outStruct.timeout = timein;
        %     outStruct.datout = datout;
        %     outStruct.statout = flagout;
            % outStruct.igap = igap;
        end
        return;
    end
end

% outStruct.igap = igap;

%% Search for the optimal order (p,q)
if isfield( inStruct, 'aka')
    % Akaike coefficient matrix given as input
    aka = inStruct.aka;
    
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
    
    % seg have ML length if mseg==0
    if mseg~=0
        if ML>=mseg
            seg = seg(1:mseg);
        end
    end
    if verbflag
        fprintf('Step 4 - Order estimation\n\nPlease wait...\n\n');
        fprintf('%d points will be used for the grid of ARMA models.\n', length(seg));
    end
    if ~outStruct.flags.auto_flag
        if exist( 'akaname', 'var' )
            if verbflag
                aka = armaord( seg, 'pmin', outStruct.params.pmin, ...
                    'pmax', outStruct.params.pmax, ...
                    'qmax', outStruct.params.qmax, 'w', akaname);
            else
                aka = armaord( seg, 'pmin', outStruct.params.pmin, ...
                    'pmax', outStruct.params.pmax, ...
                    'qmax', outStruct.params.qmax, ...
                    'verbose', false, 'w', akaname);
            end
            % Reduce aka matrix when it is larger than demanded
            ss = size(aka);
            if ss(1)>outStruct.params.pmax | ss(2)>outStruct.params.qmax
                aka = aka(1:(pmax-pmin+1), 1:(qmax+1));
            end
        % Note that, armaord requires the flag 'w' is the last one used
        elseif outStruct.flags.temp
            if verbflag
                aka = armaord( seg, 'pmin', outStruct.params.pmin, ...
                    'pmax', outStruct.params.pmax, ...
                    'qmax', outStruct.params.qmax, 'w' );
            else
                aka = armaord( seg, 'pmin', outStruct.params.pmin, ...
                    'pmax', outStruct.params.pmax, ...
                    'qmax', outStruct.params.qmax, 'verbose', false, 'w' );
            end
        else
            if verbflag
                aka = armaord( seg, 'pmin', outStruct.params.pmin, ...
                    'pmax', outStruct.params.pmax, ...
                    'qmax', outStruct.params.qmax);
            else
                aka = armaord( seg, 'pmin', outStruct.params.pmin, ...
                    'pmax', outStruct.params.pmax, ...
                    'qmax', outStruct.params.qmax, 'verbose', false);
            end
        end
        
    else
        if exist('akaname', 'var')
            if verbflag
                aka = autoarmaord( seg, 'w', akaname, 'rep', rep_lim, ...
                    'mseg', outStruct.params.mseg, ...
                    'mem', outStruct.params.mem);
            else
                aka = autoarmaord( seg, 'verbose', false, 'w', akaname, ...
                    'rep', rep_lim, 'mseg', outStruct.params.mseg, ...
                    'mem', outStruct.params.mem);
            end
        elseif temp
            if verbflag
                aka = autoarmaord( seg, 'w', 'rep', rep_lim, ...
                    'mseg', outStruct.params.mseg, ...
                    'mem', outStruct.params.mem);
            else
                aka = autoarmaord( seg, 'verbose', false, 'w', ...
                    'rep', rep_lim, 'mseg', outStruct.params.mseg, ...
                    'mem', outStruct.params.mem);
            end
        else
            if verbflag
                aka = autoarmaord( seg, 'rep', rep_lim, ...
                    'mseg', outStruct.params.mseg, ...
                    'mem', outStruct.params.mem);
            else
                aka = autoarmaord( seg, 'verbose', false, ...
                    'rep', rep_lim, 'mseg', outStruct.params.mseg, ...
                    'mem', outStruct.params.mem);
            end
        end
    end
end

outStruct.aka = aka;

% The optimal (p,q) pair is found. In case of coincidence the lower p is
% the preference
[cp, cq] = find( aka == min( min( aka ) ) );
q = cq - 1;
p = cp + outStruct.params.pmin - 1;
fprintf('\nOptimal order: [%d %d]\n\n', p, q); 

pred_lim = 4000;
faclim = floor( pred_lim/(p+q) );

if outStruct.params.facmax>faclim
    fprintf(2,'\nWarning: facmax greater than %d might produce issues\n\n', faclim);
end

%% 1st ARMA filling run

j = 1; % iteration-number

l0 = length( igap );

% flagout = flagin;
datout_tmp = datout;

% Number of gaps
numgap = l0/2;

% Merge flag is internal. It is activated automatically when necessary
merge_flag = 0;

if numgap==1
    fprintf('**Starting the gap-filling iterative process**\n\n');
    fprintf('Number of gaps: %d\n', numgap);
    
    [datout, flagout] = af_simp( outStruct, [1 1]); % [1 1] is necessary when 
    % the debug flag is activated to save file deb001.csv in folder 001
    outStruct.datout = datout;
    outStruct.statout = flagout;
    
    igap = indgap(flagout);
    l1 = length( igap );
    numgap = l1/2;
    fprintf('\nNumber of gaps remaining: %d\n', numgap);
    
else
    fprintf('**Starting the gap-filling iterative process**\n\n');
    fprintf('Total number of gaps: %d\n\n', numgap);
    while numgap>1
            
        [datout, flagout, ftc] = af_simp( outStruct, [j 1]); % [j 1] is used
        % when debug mode is activated to save file deb00j.csv in folder 
        % 001 where folder number refers to the ARMA filling section
        outStruct.datout = datout;
        outStruct.statout = flagout;
        
        % Activate the FT correction with ftc flag from af_simp
        if ftc     
            outStruct.flags.ft_corr = ftc;      
        end
        
        igap = indgap( flagout );
        outStruct.igap = igap;
        
        if isempty( igap )
            numgap = 0;
            fprintf('\nNumber of gaps: 0\n');
            break;
        end
                   
        % Number of gaps
        l1 = length( igap );
        numgap = l1/2;
        fprintf('\nNumber of gaps: %d\n\n', numgap);
        
        j = j + 1;

        % Termination condition: the number of gaps is not repeated 
        % more than twice in consecutive iterations
        if l1==l0
            if merge_flag==0
                merge_flag = 1;
            else
                merge_flag = 0;
                break;
            end
        else
            merge_flag = 0;
            l0 = l1;
        end
                
        % If the number of gaps is still greater than 1 it will merge some
        % of them and repeat the main loop
        if merge_flag
            j = 1;
            if numgap>1
    %             flagout( flagout~=1 ) = 0;
                fprintf( '\n**Reinicialization with gap merging**\n' );
                numgap0 = numgap;
                [flagout, go] = gapmerge( flagout, igap, outStruct.params.facint );
                outStruct.statout = flagout;
                if go==true
    %                 flagin( flagout==-1 ) = -1;
                    igap = indgap( flagout );
                    outStruct.igap = igap;
                    l0 = length( igap );
                    numgap = l0/2;
                    fprintf('\nMerged gaps: %d\n', numgap0-numgap );
                    fprintf('\nNumber of gaps: %d\n\n', numgap);
                else
                    fprintf(2, warning_m3);
                    break;
                end
            end
        end
    end
end

% Recover data segments that were taken out with sing
% datout(flagin==-1) = datin(flagin==-1);
flagout( flagin~=1 ) = 0;
outStruct.statout = flagout;
igap = indgap( flagout );
outStruct.igap = igap;

% if (exist('Llin','var'))
%     Llin = Llin + length(find(flagout~=0));
% end

%% 2nd ARMA filling run (the optimal order condition is relaxed)

if numgap > 0
    fprintf( '\n**Reducing ARMA(p,q) order for the remaining gaps**\n\n' );

    if numgap==1
        [datout, flagout] = af_simp( outStruct, [1 2], 'lastr_aka', true);
        % [1 2] is used when debug mode is activated to save file 
        % deb001.csv in folder 002 where folder number refers to the ARMA 
        % filling section
        outStruct.datout = datout;
        outStruct.statout = flagout;
        
        igap = indgap( flagout );
        l1 = length( igap );
        numgap = l1/2;
        fprintf('\nNumber of gaps remaining: %d\n', numgap);
        
    else
        while numgap>=1
            [datout, flagout, ftc] = af_simp(outStruct, [j 2], ...
                'lastr_aka', true);    % [j 2] is used when debug mode is 
            % activated to save file deb00j.csv in folder 002 where folder 
            % number refers to the ARMA filling section
            outStruct.datout = datout;
            outStruct.statout = flagout;
            
            % Activate the FT correction with ftc flag from af_simp
            if ftc
                outStruct.flags.ft_corr = ftc;      
            end
            
            igap = indgap(flagout);
            outStruct.igap = igap;
            
            if isempty( igap )
                numgap = 0;
                fprintf('\nNumber of gaps: 0\n');
                break;
            end
                       
            % Number of gaps
            l1 = length(igap);
            numgap = l1/2;
            fprintf('\nNumber of gaps: %d\n\n', numgap);
            
            j = j + 1;
            
            % Termination condition: the number of gaps is not repeated 
            % more than twice in cornsecutive iterations
            if l1==l0
                if merge_flag==0
                    merge_flag = 1;
                else
                    merge_flag = 0;
                    break;
                end
            else
                l0 = l1;
                merge_flag = 0;
            end
    
            % if mod(j,2)==1 && j>1
            %     datout = flipud( datout );
            %     flagout = fliplr( flagout );
            %     igap = L - igap + 1;
            %     igap = fliplr( igap );
            % end
            
            % If the number of gaps is still greater than 1 it will merge some
            % of them and repeat the main loop

            if merge_flag
                j = 1;
                if numgap>1
        %             flagout( flagout~=1 ) = 0;
                    fprintf( '\n**Reinicialization with gap merging**\n' );
                    numgap0 = numgap;
                    [flagout, go] = gapmerge( flagout, igap, outStruct.params.facint );
                    outStruct.statout = flagout;
        
                    if go==true
        %                 flagin( flagout==-1 ) = -1;
                        igap = indgap( flagout );
                        outStruct.igap = igap;
                        l0 = length( igap );
                        numgap = l0/2;
                        fprintf('\n Merged gaps: %d\n', numgap0-numgap );
                        fprintf('\nNumber of gaps remaining: %d\n\n', numgap);
                    else
                        fprintf(2, warning_m3);
                        break;
                    end
                end
            end
        end
    end
    outStruct.igap = igap;
end

%% 3rd ARMA filling run (One-sided extrap is activated if always_int is on)
% Fill gaps left previously due to any issue in armaint that set the flag 
% <go> to False.

if (outStruct.flags.always_int && numgap > 0)
    fprintf('\n**Restarting the gap-filling with one-sided extrap**\n\n');
    igap = indgap(flagout);
    outStruct.igap = igap;
    j = 1;

    if numgap==1
        [datout, flagout] = af_simp(outStruct, [1 3], '1s'); % [1 3] is used 
        % when debug mode is activated to save file deb001.csv in 
        % folder 003 where folder number refers to the ARMA filling section
        outStruct.datout = datout;
        outStruct.statout = flagout;

        % Activate the FT correction with ftc flag from af_simp
        if ftc     
            outStruct.flags.ft_corr = ftc;
        end
        
        igap = indgap( flagout );
        l1 = length( igap );
        numgap = l1/2;
        fprintf('\nNumber of gaps remaining: %d\n', numgap);
    else
        while numgap>=1
            [datout, flagout, ftc] = af_simp(outStruct, [j 3], '1s');
            % [j 3] is used when debug mode is activated to save file 
            % deb00j.csv in folder 003 where folder number refers to the 
            % ARMA filling section
            outStruct.datout = datout;
            outStruct.statout = flagout;

            % Activate the FT correction with ftc flag from af_simp
            if ftc     
                outStruct.flags.ft_corr = ftc;      
            end

            igap = indgap(flagout);
            outStruct.igap = igap;
           
            if isempty(igap)
                fprintf('\nNumber of gaps: 0\n');
                break;
            end

            % Number of gaps
            l1 = length(igap);
            numgap = l1/2;
            fprintf('\nNumber of gaps: %d\n\n', numgap);
            
            j = j + 1;
            
            % Termination condition: the number of gaps is not repeated 
            % more than twice in consecutive iterations
            if l1==l0
                if merge_flag==0
                    merge_flag = 1;
                else
                    merge_flag = 0;
                    break;
                end
            else
                l0 = l1;
                merge_flag = 0;
            end

            % if mod(j,2)==1 && j>1
            %     datout = flipud( datout );
            %     flagout = fliplr( flagout );
            %     igap = L - igap + 1;
            %     igap = fliplr( igap );
            % end

            % If the number of gaps is still greater than 1 it will merge some
            % of them and repeat the main loop

            if merge_flag
                j = 1; 
                if numgap>1
                    fprintf( '\n *Reinicialization with gap merging*\n' );
                    numgap0 = numgap;
                    [flagout, go] = gapmerge( flagout, igap, outStruct.params.facint );
                    outStruct.statout = flagout;

                    if go==true
                        igap = indgap( flagout );
                        outStruct.igap = igap;
                        l0 = length( igap );
                        numgap = l0/2;
                        fprintf('\n Merged gaps: %d\n', numgap0-numgap );
                        fprintf('\nNumber of gaps remaining: %d\n\n', numgap);
                    else
                        fprintf(2, warning_m3);
                        break;
                    end
                end
            end
        end
    end
    outStruct.igap = igap;
end

%% 4th ARMA filling run (One-sided extrap + relaxed optimal condition)

if (outStruct.flags.always_int && numgap > 0)
    fprintf('\n One-sided extrap + relaxed optimal condition \n');
    fprintf( ' Reducing ARMA order *\n\n' );
    igap = indgap(flagout);
    outStruct.igap = igap;

    if numgap==1
        [datout, flagout, ftc] = af_simp( outStruct, [1 4], ...
            'lastr_aka', true, '1s' ); % [1 4] input is used 
        % when debug mode is activated to save file deb001.csv in 
        % folder 004 where folder number refers to the ARMA filling section
        outStruct.datout = datout;
        outStruct.statout = flagout;

        % Activate the FT correction with ftc flag from af_simp
        if ftc     
            outStruct.flags.ft_corr = ftc;      
        end

        igap = indgap( flagout );
        l1 = length( igap );
        numgap = l1/2;
        fprintf('\nNumber of gaps remaining: %d\n', numgap);

        if numgap>0
            fprintf(2, warning_m3);
            fprintf(2, warning_m1);
            outStruct.flags.ft_corr = false;
        end

    else
        while numgap>=1
            [datout, flagout, ftc] = af_simp( outStruct, [j 4], ...
                'lastr_aka', true, '1s' ); % [j 4] input is used 
            % when debug mode is activated to save file deb00j.csv in 
            % folder 004

            % Activate the FT correction with ftc flag from af_simp
            if ftc
                outStruct.flags.ft_corr = ftc;      
            end

            igap = indgap( flagout );
            outStruct.igap = igap;

            if isempty( igap )
                % numgap = 0;
                fprintf('\nNumber of gaps remaining: 0\n');
                break;
            end

            % Number of gaps
            l1 = length(igap);
            numgap = l1/2;
            fprintf('\nNumber of gaps remaining: %d\n\n', numgap);

            j = j + 1;

            % Termination condition: the number of gaps is not repeated 
            % more than twice in consecutive iterations
            if l1==l0
                if merge_flag==0
                    merge_flag = 1;
                else
                    break;
                end
            else
                merge_flag = 0;
                l0 = l1;
            end

            % if mod(j,2)==1 && j>1
            %     datout = flipud( datout );
            %     flagout = fliplr( flagout );
            %     igap = L - igap + 1;
            %     igap = fliplr( igap );
            % end
            
            % If the number of gaps is still greater than 1 it will merge some
            % of them and repeat the main loop
            if merge_flag
                j = 1;
                if numgap>1
        %             flagout( flagout~=1 ) = 0;
                    fprintf( '\n *Reinicialization with gap merging*\n' );
                    numgap0 = numgap;
                    [flagout, go] = gapmerge( flagout, igap, outStruct.params.facint );
                    outStruct.statout = flagout;

                    if go==true
        %                 flagin( flagout==-1 ) = -1;
                        igap = indgap( flagout );
                        outStruct.igap = igap;
                        l0 = length( igap );
                        numgap = l0/2;
                        fprintf('\n Merged gaps: %d\n', numgap0-numgap );
                        fprintf('\nNumber of gaps remaining: %d\n\n', numgap);
                    else
                        fprintf(2, warning_m3);
                        fprintf(2, warning_m1);
                        outStruct.flags.ft_corr = false;
                        break;
                    end
                else
                    if numgap==1
                        fprintf(2, warnign_m3);
                        fprintf(2, warning_m1);
                        outStruct.flags.ft_corr = false;
                    end
                    break;
                end
            end
        end
    end
    outStruct.igap = igap;
end

%% FT correction of the ARMA interpolation
if outStruct.flags.ft_corr
    %1st iteration
    datout_corr = ftcorr(datout, flagin, 'cutoff', outStruct.params.cutoff_level);
    %2nd iteration
    datout_corr = ftcorr(datout_corr, flagin, 'cutoff', outStruct.params.cutoff_level);
else
    if isfield( outStruct.params, 'ft_corr' )
        if outStruct.params.ft_corr
            fprintf(2,'No FT correction can be applied due to the remaining gaps\n');
        end
    end
end

% Recover original data that was excluded during the interpolation
if outStruct.flags.reco_flag
    datout(flaglin==0) = datout_tmp(flaglin==0);
    if outStruct.flags.ft_corr
        datout_corr(flaglin==0) = datout_tmp(flaglin==0);
    end
end

%% Save output

if outStruct.flags.ft_corr
    outStruct.datout = datout_corr;
else
    outStruct.datout = datout;
end

outStruct.timeout = timein;
outStruct.aka = aka;
outStruct.igap = igap;
% These are not the original flags but processed ones after lincorr, sing, ...
outStruct.statin = flagin;
% These are the output flags
outStruct.statout = flagout;
outStruct.ord = [p q];
if exist("seg","var")
    outStruct.segord = seg;
end

if outStruct.flags.ascii_struct
    saveout(outStruct);
    fprintf('\n  Interpolation finished successfully.  \n');
end

cd ..

end
% END
