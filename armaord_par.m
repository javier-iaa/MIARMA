function varargout = armaord_par(S,varargin)
% Function [p,q,varargout] = armaord(S,varargin) estimates the optimal pair
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

% Version: 0.1
% Changes: Parallelized version of armaord 1.0.2.6
% Author: Javier Pascual-Granado
% $Date: 18/08/2015$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning off all

%% Preparing input data
N = length(S);

% Number of iterations for the algorithm
nit = 10; 

% nin = max(nargin,1)-1;
% nout = max(nargout,1) - 2;

% Change row into column
S = reshape(S,N,1);

% Standardization
S = (S-mean(S))./std(S);

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
r = pmax-pmin+1;

akamat = NaN(r,pmax+1);

% modo = false;
% 
% iw = find(strcmp(varargin,'w'),1);
% if ~isempty(iw),
%     nom = varargin{iw+1};
%     nomfich = strcat(nom,'.aka');   
%     
%     % Open file to recover akaike coefficients
%     if exist(nomfich,'file')==2,
%         akamat0 = dlmread(nomfich);
%         akamat0 = akamat0';
%         [a,b] = size(akamat0);
% 
%         if length(akamat(:)) > length(akamat0(:)),
%             akamat(1:a,1:b) = akamat0;
%             modo = true;
%             fichw = fopen(nomfich,'w');
%         else
%             akamat = akamat0(1:r,1:(pmax+1));
%         end
%     else
%         modo = true;
%         fichw = fopen(nomfich,'w');
%     end
% 
% end

%% Main loop
% This is the main loop where the Akaike coefficients are calculated for
% each pair (p,q)

% h = waitbar(0,sprintf(...
% 'Step 2 - Order estimation\n     Please wait...'));

parfor i=0:pmax,   
    akam = akamat(:,1+i);
        for j=pmin:pmax,
            try               
                model = armax(double(single(seg)),[j i], ...
                    'SearchMethod','lm','MaxIter',nit,'LimitError',1);
            catch
                continue;
            end            
            aka = aic(model);           
            akam(j-pmin+1) = aka;            
        end
        akamat(:,1+i) = akam;       
end

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
