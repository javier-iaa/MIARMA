function varargout = armaord(S, varargin)
% Function varargout = armaord(S,varargin) estimates the optimal pair
% order (p,q) for an ARMA model fitting the time series S.
% These calls are possible:
%   armaord(S)
%   armaord(S, 'w')
%   armaord(S, 'w', filename) read/write akaike coefficient matrix in filename
%   armaord(S, 'pmin', pmin) use pmin as min value for iterations
%   armaord(S, 'pmax', pmax) use pmax as max value for iterations
%   armaord(S, 'pmin', pmin, 'pmax', pmax, 'w', filename)
%
%   Outputs:
%   varargout{1} = Akaike matrix
%   varargout{2} = pmin
%   varargout{3} = pmax
%   varargout{4} = size of the segment evaluated
%
% By Javier Pascual-Granado
% <a href="matlab:web http://www.iaa.es;">IAA-CSIC, Spain</a>
%
% Version: 1.2.1
%
% Changes:
% - Default values for pmin, pmax, qmax have changed.
% - Fixed a bug.
%
% Date: 9/3/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Information Criterion
IC = 'AICc';

warning off all

%% Preparing input data

% Control verbose: 1 (default) means verbose mode, 0 is silent mode
iver = find(strcmp(varargin,'verbose'), 1);
if isempty(iver),
    verbflag = 1;
else
    verbflag = varargin{iver+1};
end

N = length(S);

switch IC

    case 'AIC'
        ext = '.aka';
        
    case 'AICc'
        ext = '.akc';
        
    case 'BIC'
        ext = '.bic';
        
    case 'FPE'
        ext = '.fpe';
        
    case 'HQ'
        ext = '.hq';
        
    otherwise
        error('Selection criterion is not defined properly')
end
    
% Change row into column
S = reshape(S, N, 1);

% Standardization
S = stnorm(S);

% Values for pmin, pmax, qmax.
ipm = find(strcmp(varargin,'pmin'), 1);
if isempty(ipm),
    pmin = 1; % Default
else
    pmin = varargin{ipm+1};
end
ipM = find(strcmp(varargin,'pmax'), 1);
if isempty(ipM),
    pmax = 20; % Default
else
    pmax = varargin{ipM+1};
end
ipq = find(strcmp(varargin,'qmax'), 1);
if isempty(ipq),
    qmax = pmax;
else
    qmax = varargin{ipq+1};
end

seg = S;

%% Akaike matrix
r = pmax - pmin + 1;
c = qmax + 1;
nmod = (pmax-pmin+1)*(qmax+1);

akamat = NaN( r, c );

fprintf('\n Coefficient p in range [%d,%d] and q in [0,%d]\n', pmin, pmax, qmax);
fprintf(' %d models to be calculated\n', nmod);

nsav = 0;

iw = find( strcmp(varargin,'w'), 1 );
if ~isempty(iw),    
    if length(varargin)>iw
        idname = varargin{iw+1};
    else
        idname = 'temp';
    end
    idname = [ idname '_' num2str(N) ];
    pminstr = num2str( pmin, '%03.f' );
    pmaxstr = num2str( pmax, '%03.f' );
    qmaxstr = num2str( qmax, '%03.f' );
    headline = [pminstr '_' pmaxstr '_' qmaxstr '\n' ]; % Header info for aka file
    nomfich = [idname ext];
    
    % Look for aka files, get the parameters and load akamat
    d = dir( nomfich );
      
    if ~isempty(d)
        fichr = fopen(nomfich,'r');
        k=1;
        l = fgetl(fichr);
        while ~isempty(l)
            f{k} = l;
            k=k+1;
            l = fgetl(fichr);
        end
        fclose(fichr);
        nf = length(f);
        pmax0 = zeros(1,nf);
        pmin0 = zeros(1,nf);
        qmax0 = zeros(1,nf);
    
        for k=1:nf
            lst_cells = regexp( f{k}, '_', 'split');
            pmin0(k) = str2num( lst_cells{end-2} );
            pmax0(k) = str2num( lst_cells{end-1} );
            qmax0(k) = str2num( lst_cells{end} );
        end
        
        % Select the file with highest pmax, then lowest pmin, the highest 
        % qmax, i.e. the file with the broadest range. Note that this might 
        % be not the best strategy in some cases so this might change in 
        % future versions of the program.
        pmx0 = max( pmax0 );
        pmn0 = min( pmin0( pmax0==pmx0 ) );
        qmx0 = max( qmax0( pmax0==pmx0 & pmin0==pmn0 ) );
        pmax0 = pmx0;
        pmin0 = pmn0;
        qmax0 = qmx0;
               
        fprintf(' Found Akaike file %s\n', nomfich);
        akamat0 = dlmread( nomfich, '', nf+1, 0 );
        akamat0 = akamat0';
              
        if ( pmin >= pmin0 )
            if ( pmin > pmax0 )
                % Case ix) 
                % Disjoint intervals. Not implemented yet
                err = MException('ResultChk:OutOfRange', ...
                    'Outside of expected range, not implemented yet');
                throw(err)
            elseif ( pmax <= pmax0 )
                if ( qmax <= qmax0 )
                    % Case iv) 
                    ii0 = pmin - pmin0 + 1;
                    ie0 = pmax - pmin0 + 1;
                    je0 = qmax + 1;
                    akamat = akamat0( ii0:ie0, 1:je0 );
                    varargout{1} = akamat;
                    fprintf(' All coeffs. already calculated\n');
                    return
                else
                    % Case iii)
                    % Coeffs. in the subset pmax x je are read from akamat0
                    ii0 = pmin - 1;
                    ie0 = pmax - 1;
                    je = qmax0 + 1;
                    akamat(:,1:je ) = akamat0( ii0:ie0, :);
                    nsav = pmax*je;
                end
            else
                if (qmax<=qmax0)
                    % Case ii)
                    % Coeffs. in the subset ie x qmax are read from akamat0
                    ie = pmax0 - pmin + 1;
                    ii0 = pmin - 1;
                    je0 = qmax + 1;
                    akamat(1:ie ,:) = akamat0( ii0:end, 1:je0);
                    nsav = ie*(qmax+1);
                else
                    % Case i)
                    % Coeffs. in subset ie x je are read from akamat0
                    ie = pmax0 - pmin + 1;
                    ii0 = pmin - pmin0 + 1;
                    je = qmax0 + 1;
                    akamat(1:ie ,1:je) = akamat0( ii0:end, :);
                    nsav = ie*je;
                end
            end
        else
            if (pmax < pmin0)
                % Case x)
                % Disjoint intervals. Not implemented yet
                err = MException('ResultChk:OutOfRange', ...
                    'Outside of expected range, not implemented yet');
                throw(err)
            elseif ( pmax < pmax0 )
                if ( qmax <= qmax0 )
                    % Case viii)
                    % Coeffs. in the subset (pmax -ii) x qmax are read from 
                    % akamat0
                    ii = pmin0 - pmin + 1;
                    ie0 = pmax - pmin0 + 1;
                    je0 = qmax + 1;
                    akamat( ii:end, : ) = akamat0( 1:ie0, 1:je0);
                    nsav = (pmax - ii)*qmax;
                else
                    % Case vii)
                    % Coeffs. in the subset (pmax -ii) x je are read from 
                    % akamat0
                    ii = pmin0 - pmin + 1;
                    ie0 = pmax - pmin0 + 1;
                    je = qmax0 + 1;
                    akamat( ii:end, 1:je ) = akamat0( 1:ie0, :);
                    nsav = (pmax - ii)*je;
                end
            else
                if (qmax<=qmax0)
                    % Case vi)
                    % Coeffs. in the subset (ie-ii+1) x qmax are read from 
                    % akamat0
                     ii = pmin0 - pmin + 1;
                     ie = pmax0 - pmin + 1;
                     je0 = qmax + 1;
                     akamat(ii:ie ,:) = akamat0( :, 1:je0);
                     nsav = (ie - ii + 1)*qmax;
                else
                    % Case v)
                    % Coeffs. in subset (ie-ii+1) x je are read from akamat0
                    ii = pmin0 - pmin + 1;
                    ie = pmax0 - pmin + 1;
                    je = qmax0 + 1;
                    akamat(ii:ie ,1:je) = akamat0;
                    nsav = (ie - ii + 1)*je;
                end
            end
        end
    end
    
    ncal = nmod - nsav;
    fprintf(' %d models loaded, %d need to be calculated.\n\n', nsav, ncal);
    fichw = fopen(nomfich, 'w');
    
    % Write header lines
     if ~isempty(d)
         for k=1:nf
             fprintf(fichw, [ f{k} '\n' ]);
         end
     end
    fprintf(fichw, [headline '\n']);
end

%% Main loop
% This is the main loop where the Akaike coefficients are calculated for
% each pair (p,q)

tee = 0;

for i=0:qmax,
    tic;    
    akam = akamat(:,1+i);
    nmod = 0;

    for j=pmin:pmax
        ii = j - pmin + 1;
        aka = akam(ii);

        if ~isnan(aka)
            if ~isempty(iw),
                fprintf( fichw, '%f ', aka);
            end
            continue;
        end

        try               
            model = armax(seg, [j i], 'SearchMethod', 'lsqnonlin');
            nmod = nmod + 1;
        catch E
            % Report the error messages to a file
            if exist('nomfich','var'),
                errfich = strcat( [ nomfich(1:end-4) '_EM.htm' ] );
                fich = fopen(errfich,'a');
            else
                fich = fopen('Error_Messages.htm','a');
            end
            msg = getReport(E);
            fprintf(fich,'<br>(%d,%d) order failed<br/>',j,i);
            fprintf(fich,'%s\n',msg);
            fclose(fich);

            if ~isempty(iw)
                fprintf(fichw,'NaN ');
            end
            continue;
        end

        aka = aicplus(model, IC);

        akam(ii) = aka;

        if ~isempty(iw)
            fprintf( fichw, '%f ', aka);
        end 
    end

    akamat(:, 1+i) = akam;

    if ~isempty(iw)
        fprintf( fichw, '\n');
    end

    per = 100*(i+1)/(1+qmax);
    ee = toc;
    tee = tee + ee;
    if verbflag
        fprintf(' %.f%%...  %.f s last iteration, %.f s since start, %d models calculated\n', ...
            per, ee, tee, nmod);
    end
end

if ~isempty(iw),
    fclose(fichw);
end

%% Optional outputs
if nargout==1,
    varargout{1} = akamat;
elseif nargout==2,
    varargout{1} = akamat;
    varargout{2} = L;
elseif nargout==3,
    varargout{1} = akamat;
    varargout{2} = pmin;
    varargout{3} = pmax;
elseif nargout==4,
    varargout{1} = akamat;
    varargout{2} = pmin;
    varargout{3} = pmax;
    varargout{4} = L;
elseif nargout==5,
    varargout{1} = akamat;
    varargout{2} = pmin;
    varargout{3} = pmax;
    varargout{4} = L;
    varargout{5} = qmax;
end
