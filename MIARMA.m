function strout = MIARMA(varargin)
% function strout = MIARMA(varargin) 
% interpolates datapoints in a gapped time series using ARMA models to
% predict the segments of data.
% Inputs:   filename of an ASCII file containing the input data array, the
%           corresponding time array, and the status array
%           MIARMA( filename )
%           Also a parameter filename can be given as second input
%           MIARMA( filename, param )
%
% Alternative call:
%           MIARMA( strdata )
%            where strdata is a struct that contains the necessary inputs:
%             time, data, stat
%            and a set of optional parameters:
%             igap - is the gap indexes array
%             aka - is the Akaike coefficient matrix
%             temp - boolean, 1 to save temp files 0 otherwise
%             params - a struct that contain constant parameters simpl,
%              temp, facmin, facmax, npi, npz, repmax, pmin, pmax, mseg, 
%              nuc, and always_int. 
%             Consult the documentation for a detailed description of each 
%             of these parameters.
%
% Outputs:  If no output variable is given an ASCII file is created 
%           containing the output data and the corresponding time and
%           status
%
%           For test purposes any of theses calls are available:
%           strout = MIARMA( filename )
%           which is the output structure. It might contains timeout,
%            datout, statout    
%           and also the following variables:
%            aka - Akaike coefficient matrix
%            ord - optimal order (p,q)
%            igap - gap indexes
%            params - a complete list of parameters used for computation
%
% By Javier Pascual-Granado
% <a href="matlab:web http://www.iaa.es;">IAA-CSIC, Spain</a>

% Dependencies:     armaord_par.m   0.1
%                   armaord.m       1.0.4
%                   indgap.m        1.0.2
%                   armafill.m      1.2.1
%                   lincorr.m       1.0.4
%                   sing.m
%                   af_simp.m       0.1.4
%                   polintre.m
%                   armaint.m       1.3.4
%                   pred.m          1.0.1
%
% Version: 1.5.9
%
% Changes: - waitbar substituted by percentages
%
% Date: 04/12/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off all

%% Some definitions
version = '1.5.9';

lgaps0 = NaN;
Llin = NaN;

% Call definition of cellfind function
cellfind = @(string)(@(cell_contents)(strcmp(string,cell_contents)));

fprintf('Warning: All gaps must be correctly flagged for the gap-filling\n')
fprintf('algorithm to give an adequate output\n\n');

%% Input data
if ischar( varargin{1} )
    
    filename = varargin{ind};
    
    % This is for CoRoT data only    
    if strcmp( filename(1:8), 'AN2_STAR' )
        field = 'sismo';
    elseif strcmp( filename(1:8), 'EN2_STAR' )
        field = 'exo';
    else
        field = 'unknown';
    end
    
    % Here data is imported from an ASCII file having 3 columns: time, flux
    % and status
    strdata = importdata(filename);
    timein = strdata.data(:,1);
    datin = strdata.data(:,2);
    flagin = strdata.data(:,3);
    
else
    strdata = varargin{1};
    timein = strdata.time;
    datin = strdata.data;
    flagin = strdata.stat;
    field = 'unknown';
    
end

datout = datin;
timeout = timein;
L = length(datin);

%% Parameters
% Optional parameter input file
if (nargin == 2 && ischar(varargin{2}))
    parname = varargin{2};
    pardata = importdata(parname);

    % Flag to activate (or not) ascii output
    pari = strcmp(pardata.textdata,'ascii');
    pari = circshift(pari,-2);
    if ~isempty(pardata.data(pari))
        ascii_struct = pardata.data(pari);
    end
    
    % Flag to decide between af_simp.m and armafill.m
    pari = strcmp(pardata.textdata,'simpl');
    pari = circshift(pari,-2);
    if ~isempty(pardata.data(pari))
        simpl = pardata.data(pari);
    end

    % Default values of the parameters.
    % Min. ratio between segment length and number of parameters for
    % the model. facmin must be >= 3
    pari = strcmp(pardata.textdata,'facmin');
    pari = circshift(pari,-2);
    if ~isempty(pardata.data(pari))
        facmin = pardata.data(pari);
    end

    % Max. ratio between segment length and number of parameters for
    % the model. facmax = 6;
    pari = strcmp(pardata.textdata,'facmax');
    pari = circshift(pari,-2);
    if ~isempty(pardata.data(pari))
        facmax = pardata.data(pari);
    end

    % Lower limit in gap length for the ARMA interpolation
    % (below this limit linear interpolation is used)
    pari = strcmp(pardata.textdata,'npi');
    pari = circshift(pari,-2);
    if ~isempty(pardata.data(pari))
        npi = pardata.data(pari);
    end

    % Lower limit in data segment length for the ARMA interpolation.
    % This must be at least d*facmin and, as min(d)=min(p+q)=pmin+0,
    % npz must be at least pmin*facmin
    pari = strcmp(pardata.textdata,'npz');
    pari = circshift(pari,-2);
    if ~isempty(pardata.data(pari))
        npz = pardata.data(pari);
    end

    % Number of iterations of the gap-filling loop
    pari = strcmp(pardata.textdata,'repmax');
    pari = circshift(pari,-2);
    if ~isempty(pardata.data(pari))
        repmax = pardata.data(pari);
    end
    
    % Maximum length of the segment used to calculate ARMA order
    pari = strcmp(pardata.textdata,'mseg');
    pari = circshift(pari,-2);
    if ~isempty(pardata.data(pari))
        mseg = pardata.data(pari);
    end
    
    % Range for the search of the optimal ARMA orders [pmin,pmax] 
    % The MA order is search in the range [0,p]
    pari = strcmp(pardata.textdata,'pmin');
    pari = circshift(pari,-2);
    if ~isempty(pardata.data(pari))
        pmin = pardata.data(pari);
    end
    pari = strcmp(pardata.textdata,'pmax');
    pari = circshift(pari,-2);
    if ~isempty(pardata.data(pari))
        pmax = pardata.data(pari);
    end

    % Number of workers to used for parallelization. Default is 1 meaning
    % no parallelization.
    % The max number of workers permitted in the standard version of Matlab 
    % is 12, therefore, nuc should be less than 12 if MIARMA is not running
    % under Matlab Distributed Computing.
    pari = strcmp(pardata.textdata,'nuc');
    pari = circshift(pari,-2);
    if ~isempty(pardata.data(pari))
        nuc = pardata.data(pari); 
    end
    
    % Always interpolate or not
    pari = strcmp(pardata.textdata,'always_int');
    pari = circshift(pari,-2);
    if ~isempty(pardata.data(pari))
        always_int = pardata.data(pari);
    end

    % Decides wheter to save the Akaike matrix at a temp file
    pari = strcmp(pardata.textdata,'temp');
    pari = circshift(pari,-2);
    if ~isempty(pardata.data(pari))
        temp = pardata.data(pari);
    end
    
else
    % default values
    simpl = true;
    temp = false;
    always_int = true;
    nuc = 1;
    mseg = 1000;
    pmin = 2;
    pmax = 30;
    repmax = 3;
    npz = 8;
    npi = 6;
    facmax = 6;
    facmin = 4;
    ascii_struct = false;

    if isfield( strdata, 'params' )
        
        if isfield( strdata.params, 'simpl')
            simpl = strdata.params.simpl;
        end

        if isfield( strdata.params, 'facmin')
            facmin = strdata.params.facmin;
        end

        if isfield( strdata.params ,'facmax')
            facmax = strdata.params.facmax;
        end
    
        if isfield( strdata.params, 'npi')
            npi = strdata.params.npi;
        end
    
        if isfield( strdata.params, 'npz' )
            npz = strdata.params.npz;
        end
    
        if isfield( strdata.params, 'repmax')
            repmax = strdata.params.repmax;
        end

        if isfield( strdata.params, 'pmin')
            pmin = strdata.params.pmin;
        end

        if isfield( strdata.params, 'pmax')
            pmax = strdata.params.pmax;
        end

        if isfield( strdata.params, 'mseg')
            mseg = strdata.params.mseg;
        end    
    
        if isfield( strdata.params, 'nuc' )
            nuc = strdata.params.nuc;
        else
            nuc = 1;
        end

        if isfield( strdata.params, 'always_int' )
            always_int = strdata.params.always_int;
        end

        if isfield( strdata.params, 'temp' )
            temp = strdata.params.temp;
        end
        
        if isfield(strdata.params, 'ascii_struct')
            ascii_struct = strdata.params.ascii_struct;
        end
        
    end
    
end

% Save parameters used in the computation for transparency
strout.params.simpl = simpl;
strout.params.temp = temp;
strout.params.always_int = always_int;
strout.params.nuc = nuc;
strout.params.mseg = mseg;
strout.params.pmin = pmin;
strout.params.pmax = pmax;
strout.params.repmax = repmax;
strout.params.npz = npz;
strout.params.npi= npi;
strout.params.facmax = facmax;
strout.params.facmin = facmin;
strout.params.ascii_struct = ascii_struct;

% List of parameters for armafill/af_simp
% params = [facmin facmax npi pmin rstd];
params = [facmin facmax npi pmin];

% Other optional parameters
if isfield( strdata, 'igap')
    % Indexes given as input
    igap = strdata.igap;
end

if isfield( strdata, 'aka')
    % Akaike coefficient matrix given as input
    aka = strdata.aka;
end

%% Building the gap indexes
% Apply cellfind to each cell of cell array varargin
if isempty(find(cellfun(cellfind('igap'),varargin),1))
    % First step: valid data is flagged zero
    % This is for CoRoT data only 
    if strcmp(field,'exo')
        flagin(flagin==8 | flagin==16 | flagin==64 | flagin==128 | ...
        flagin==24 | flagin==72 | flagin==136 | flagin==80 | ...
        flagin==144 | flagin==192 | flagin==88 | flagin==200 | ...
        flagin==152 | flagin==208) = 0;    
    elseif strcmp(field,'sismo')
        flagin(flagin==64 | flagin==256 | flagin==512 | flagin==320 | ...
        flagin==578) = 0;
    end
        
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
%         [datout,flagin] = lincorr(datin,flagin,igap,npi,rstd);
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
%             L = length(datout);
            Llin = length(find(flagout~=0));
            if exist('filename','var')
                fout = cat(2,filename(1:end-4),'.agfs');
            else
                fout = 'curve.agfs';
            end
            fich = fopen(fout,'w');
    
            fprintf(fich,'# code: MIARMA\n');
            fprintf(fich,'# version: %s\n',version);
            fprintf(fich,'# model (p,q): not calculated\n');
            fprintf(fich,'# length: %d\n',L);
            fprintf(fich,'# gaps_arma: 0\n');
            fprintf(fich,'# gaps_linear: %d\n',Llin);
            fprintf(fich,'# facmin: %d\n',facmin);
            fprintf(fich,'# facmax: %d\n',facmax);
            fprintf(fich,'# npi: %d\n',npi);
            fprintf(fich,'# npz: %d\n',npz);
            fprintf(fich,'# niter: %d\n',repmax);
            fprintf(fich,'# mseg: %d\n',mseg);
            fprintf(fich,'x y z\n');
   
            for i=1:L
                fprintf(fich,'%16.12f %16.13f %d\n',...
                    timeout(i),datout(i),flagout(i));
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
if isempty(find(cellfun(cellfind('aka'),varargin),1))
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
    
    % Mseg is the max. length of the segment use to calculate the order
    if ML>mseg
        seg = seg(1:mseg);
    end
    
    fprintf('Step 4 - Order estimation\n    Please wait...\n');
    
    % Optimal order (p,q) for the ARMA model of seg
    if nuc == 1
        if exist( 'filename', 'var' )
            fileaka = [filename(1:end-4) '_' num2str(pmin) ...
                '_' num2str(pmax)];
            aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, 'w', fileaka);
            
        elseif temp            
            fileaka = ['temp_' num2str(pmin) '_' num2str(pmax)];
            aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, 'w', fileaka);
            
        else            
            aka = armaord( seg, 'pmin', pmin, 'pmax', pmax);            
        end
        
    else

        % Open a pool for parallel computation with nuc workers.
        
        % Version control
        rel = version;
        rnuml = rel(end-5:end-1);
        rnum = rnuml(1:end-1);
        num = str2double( rnum );
        
        if num > 2013 || strcmp(rnuml, '2013b')
            parpool(nuc);
        else
            matlabpool(nuc);
        end

    %     if exist('filename','var'),
    %         aka = armaord_par(seg,'pmin',pmin,'pmax',pmax,'w',filename(1:end-4));
    %     else
             aka = armaord_par(seg,'pmin',pmin,'pmax',pmax);
    %     end
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
l1 = l0-1;

flagout = flagin;

% Number of gaps
numgap = floor(l0/2);

if numgap==1
    if simpl,
        [datout, flagout] = af_simp( datout, flagout, aka, igap, params, j );
    else
        [datout, flagout] = armafill( datout, flagout, aka, igap, params, j );
    end
    
    if mod(j,2)==0
        datout = datout(end:-1:1);
    end
    
    % Recover data segments that were taken out with sing
%     datout(flagin==-1) = datin(flagin==-1);

%     datout(flagout~=0) = [];
%     timeout(flagout~=0) = [];
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
                l1 = 0;
                break;
            end
            
            flagout = sing(flagout, npz, igap);

            igap = indgap(flagout,1);
            
            l1 = length(igap);

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

        numgap = floor(l1/2);

        if numgap>1
            [flagout,go] = gapmerge(flagout,igap);
            if go==true
                flagin(flagout==-1)=-1;
                igap = indgap(flagout,j+1);
                rep = 0;
                l0 = length(igap);
                l1 = l0-1;
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

if ~ascii_struct
    % This occurs only when MIARMA is called for test purposes only
    strout.timeout = timeout;
    strout.datout = datout;
    strout.statout = flagout;
else
    % This is the situation in a standard call
    if exist('filename','var')
        fout = cat(2,filename(1:end-4),'.agfs');
    else
        fout = 'curve.agfs';
    end
    fich = fopen(fout,'w');
    
    fprintf(fich,'# code: MIARMA\n');
    fprintf(fich,'# version: %s\n',version);
    fprintf(fich,'# model (p,q): %d %d\n',p,q);
    fprintf(fich,'# length: %d\n',L);
    fprintf(fich,'# gaps_arma: %d\n',lgaps0-Llin);
    fprintf(fich,'# gaps_linear: %d\n',Llin);
    fprintf(fich,'# facmin: %d\n',facmin);
    fprintf(fich,'# facmax: %d\n',facmax);
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