function varargout = armaord(S, varargin)
% Function varargout = armaord(S,varargin) estimates the optimal pair
% order (p,q) for an ARMA model fitting the time series S.
% These calls are possible:
%   armaord(S)
%   armaord(S,'w',filename) read/write akaike coeffient matrix in filename
%   armaord(S,'pmin',pmin) use pmin as min value for iterations
%   armaord(S,'pmax',pmax) use pmax as max value for iterations
%   armaord(S,'pmin',pmin,'pmax',pmax,'w',filename)
%
%   Outputs:
%   varargout{1} = Akaike matrix
%   varargout{2} = pmin
%   varargout{3} = pmax
%   varargout{4} = size of the segment evaluated

% Version: 1.0.10
% Changes: bug fixes.
% Author: Javier Pascual-Granado
% $Date: 21/05/2020$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Information Criterion
IC = 'AICc';

warning off all

%% Preparing input data
N = length(S);

% Number of iterations for the algorithm
% nit = 10; 

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
    

% nin = max(nargin,1)-1;
% nout = max(nargout,1) - 2;

% Change row into column
S = reshape(S, N, 1);

% Standardization
S = stnorm(S);

% If the parameters pmin, pmax are not introduced the subroutine uses
% an estimate based on the length of the data segment
ipm = find(strcmp(varargin,'pmin'), 1);
if isempty(ipm),
    pmin = round(N/3-1);
else
    pmin = varargin{ipm+1};
end

ipM = find(strcmp(varargin,'pmax'), 1);
if isempty(ipM),
    pmax = round(N/2-1);
else
    pmax = varargin{ipM+1};
end

seg = S;

%% Akaike matrix
r = pmax - pmin + 1;

akamat = NaN( r, pmax+1 );

% We use modo as a flag to indicate whether to write the results on a file
% or not
modo = false;

iw = find( strcmp(varargin,'w'), 1 );
if ~isempty(iw),
    nom = varargin{iw+1};
    nomfich = strcat( nom, ext );
    
    % Look for aka temp files, get the p interval and load akamat
    idname = nom(1:end-8);
    d = dir([idname '*' ext]);
%     d = dir(['temp_*_*' ext]);
    pmax0 = zeros(size(d));
    pmin0 = zeros(size(d));
    if ~isempty(d),
        for k=1:length(d)
            lst_cells = regexp( d(k).name, '_', 'split');
            pmin0(k) = str2num( lst_cells{end-1} );
            cell_pmax= regexp( lst_cells{end}, '\.', 'split');
            pmax0(k) = str2num( cell_pmax{1} );
        end
        pmx0 = max(pmax0);
        pmin0 = pmin0(pmax0==pmx0);
        pmax0 = pmx0;
%         lst_cells = regexp( d(end).name, '_', 'split');
%         pmin0 = str2num( lst_cells{2} );
%         cell_pmax = lst_cells{3};
%         pmax0 = str2num( cell_pmax(1:end-4) );
%         akamat0 = dlmread( d(end).name );
        pmin0str = num2str(pmin0,'%03.f');
        pmax0str = num2str(pmax0,'%03.f');
        akaname = [idname '_' pmin0str '_' pmax0str ext];
        akamat0 = dlmread( akaname );
        akamat0 = akamat0';
        [a, b] = size( akamat0 );
    
        if (pmin >= pmin0)
            if (pmin > pmax0)
                % Disjoint intervals. Not implemented yet
                err = MException('ResultChk:OutOfRange', ...
                    'Outside of expected range, not implemented yet');
                throw(err)
            elseif (pmax <= pmax0)
                akamat = akamat0(1:r,1:(pmax+1));
                varargout{1} = akamat;
                return
            else
                % Note that the final Akaike matrix collect the maximum
                % number of coefficients so in this case it will be a
                % pmin0 x pmax matrix
                pmaxstr = num2str(pmax,'%03.f');
                nomfich = [idname '_' pmin0str '_' pmaxstr ext];
                r = pmax - pmin0 + 1;
                akamat = NaN( r, pmax+1 );
                akamat(1:a,1:b) = akamat0;
            end
        else
            if (pmax < pmin0)
                % Disjoint intervals. Not implemented yet
                err = MException('ResultChk:OutOfRange', ...
                    'Outside of expected range, not implemented yet');
                throw(err)
            elseif (pmax <= pmax0)
                pminstr = num2str(pmin,'%03.f');
                nomfich = [idname '_' pminstr '_' pmax0str ext];
                r = pmax0 - pmin + 1;
                akamat = NaN( r, pmax0+1 );
            end
            % if pmax > pmax0 nothing needs to be done, nomfich and akamat 
            % are preserved
        end
    end
    modo = true;
    fichw = fopen(nomfich, 'w');
end

%% Main loop
% This is the main loop where the Akaike coefficients are calculated for
% each pair (p,q)

% h = waitbar(0,sprintf(...
% 'Step 2 - Order estimation\n     Please wait...'));

for i=0:pmax,
    
%     waitbar((i+1)/(1+pmax),h)
    fprintf('%.f%%... ', 100*(i+1)/(1+pmax));
    
    akam = akamat(:,1+i);
    
%     ir = i*r;
    
%     fprintf('Step %d\n',1+i);
        for j=pmin:pmax,
            ii = j-pmin+1;
%             jj = 1+i;
%             ind = ir + ii;
%             aka = akamat(ii,jj);
            aka = akam(ii);
%             aka = akamat(ind);
            if ~isnan(aka),
                if modo,
                    fprintf(fichw,'%f ',aka);
                end
%                 fprintf(repmat('\b',1,5));
                continue;
            end
            
%             d = i+j;     

            % Special four-stage LS-IV iterative algorithm
%             seg1 = seg;
            try               
%                 model = armax(double(single(seg)),[j i],'Focus', ...
%                     'stability','SearchMethod','lm','MaxIter',nit, ...
%                     'LimitError',1);
%                 model = armax(seg, [j i]);
                model = armax(seg, [j i], 'SearchMethod', 'lsqnonlin');
            catch E
                % Report the error messages to a file
                if exist('nomfich','var'),
                    errfich = strcat(nom,'_EM.htm');
                    fich = fopen(errfich,'a');
                else
                    fich = fopen('Error_Messages.htm','a');
                end
                msg = getReport(E);
                fprintf(fich,'<br>(%d,%d) order failed<br/>',j,i);
                fprintf(fich,'%s\n',msg);
                fclose(fich);
                
                if modo,
                    fprintf(fichw,'NaN ');
                end
%                 fprintf(repmat('\b',1,5));
                continue;
            end
            
            aka = aicplus(model,IC);
            
            akam(ii) = aka;
%             akamat(ii,jj) = aka;
            
            if modo,
                fprintf(fichw,'%f ',aka);
            end 
        end
        
        akamat(:,1+i) = akam;
        
        if modo,
            fprintf(fichw,'\n');
        end
%         fprintf(repmat('\b',1,5));
end

% close(h);

if modo,
    fclose(fichw);
end

% [cp,cq]=find(akamat==min(min(akamat)));
% 
% % The optimal (p,q) pair is found. In case of coincidence the lower p is
% % the preference
% q = cq - 1;
% p = cp + pmin - 1;

%% Optional outputs
% This part of the code is just for test purposes
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
end
