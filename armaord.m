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
%   varargout{2} = ppmax
%   varargout{3} = qpmax
%   varargout{4} = size of the segment evaluated
%
% where ppmax and qqmax are the pmax and qmax re-evaluated from the akaike
% coefficient matrix in filename
%
% By Javier Pascual-Granado
% <a href="matlab:web http://www.iaa.es;">IAA-CSIC, Spain</a>
%
% Version: 2.2
%
% Changes:
% - Parfor loop is now calculating each model sequentially instead of row
% by row. This allows more parallelization of the loop and increase the
% performance.
% - Other minor fixes.
%
% Date: 14/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Default values and initial setup

% Information Criterion (IC)
% Change IC value to: 'AIC', 'BIC', 'FPE', 'HQ' for other criteria
% If necessary, change the extension to: '.aka', '.bic', '.fpe', '.hq'
IC = 'AICc';
ext = '.akc';

warning('off', 'all');
verbflag = 1;
pmin = 1;
pmax = 20;
qmax = 20;

% Control verbose: 1 (default) means verbose mode, 0 is silent mode
iver = find(strcmp(varargin,'verbose'), 1);
if ~isempty(iver)
    verbflag = varargin{iver+1};
end

% Preprocess input data
S = reshape(S, [], 1);

% Standardization
S = stnorm(S);
N = length(S);

% Override default pmin, pmax, qmax if specified
ipm = find(strcmp(varargin,'pmin'), 1);
if ~isempty(ipm)
    pmin = varargin{ipm+1};
end
ipM = find(strcmp(varargin,'pmax'), 1);
if ~isempty(ipM)
    pmax = varargin{ipM+1};
end
ipq = find(strcmp(varargin,'qmax'), 1);
if ~isempty(ipq)
    qmax = varargin{ipq+1};
end

%% Akaike matrix preallocation
r = pmax - pmin + 1;
c = qmax + 1;
akamat = NaN(r, c);

fprintf('Coefficient p in range [%d,%d] and q in [0,%d]\n', pmin, pmax, qmax);

% Optional file handling and load previous computations
iw = find(strcmp(varargin,'w'), 1);

models_to_calc = r * c;

headline = sprintf('%03.f_%03.f_%03.f\n',pmin, pmax, qmax);

if ~isempty(iw)
    idname = 'temp';
    if length(varargin) > iw
        idname = varargin{iw+1};
    end

    nomfich = sprintf('%s_%d%s', idname, N, ext);    
    if isfile(nomfich)
        fprintf('\nLoading previous calculations from %s\n', nomfich);
        
        % Read the header
        fid = fopen(nomfich, 'r');
        k=1;
        f = fgetl(fid);
        while ~isempty(f)
            header{k} = f; %#ok<AGROW>
            k = k+1;
            f = fgetl(fid);
        end
        nh = length(header);
        
        % Read the data (skip the header)
        existingData = readmatrix(nomfich, 'NumHeaderLines',nh+1,'FileType','text');
        
        % Load existing data into akamat
        loaded_rows = size(existingData, 1);
        loaded_cols = size(existingData, 2);
        akamat(1:loaded_rows, 1:loaded_cols) = existingData;
        
        % Count loaded models
        models_loaded = sum(~isnan(existingData(:)));
        models_to_calc = models_to_calc - models_loaded;

        if models_to_calc<=0
            fprintf(' %d models loaded, 0 need to be calculated.\n\n', models_loaded);
            ppmax = loaded_rows + pmin - 1;
            qpmax = loaded_cols - 1;
            varargout = {akamat, ppmax, qpmax, N};
            fprintf('Extended Akaike matrix with pmax %d and qmax %d\n', ppmax, qpmax);
            return
        else
            fprintf(' %d models loaded, %d need to be calculated.\n\n', models_loaded, models_to_calc);
        end
 
        fichw = fopen(nomfich, 'w');

        for k=1:nh
            fprintf(fichw, [ header{k} '\n']);
        end

        % Extract model indexes from last row
        C = textscan(header{end},'%03.f_%03.f_%03.f');
        if C{2}<pmax | C{3}<qmax
            fprintf(fichw, headline);
        end
        fprintf(fichw, '\n');

    else
        fid = fopen(nomfich, 'w');

        % Header info for aka file
        fprintf(fid, headline);
        fprintf(fid, '\n');
        fclose(fid);
    end
end

%% Main parallel loop to calculate missing Akaike coefficients

tic;
% Identify all (i,j) combinations that still need to be calculated
[I, J] = ndgrid(0:qmax, pmin:pmax);

% Convert (i,j) to the corresponding position in akamat
row = J - pmin + 1;
col = I + 1;

mask = isnan(akamat(sub2ind(size(akamat), row, col)));

% List of pending combinations
pending_i = I(mask);
pending_j = J(mask);

nPending = numel(pending_i);

% Calculate only the missing models
values = NaN(nPending, 1);

parfor k = 1:nPending
    i = pending_i(k);
    j = pending_j(k);

    try
        model = armax(S, [j i], 'SearchMethod', 'lsqnonlin');
        values(k) = aicplus(model, IC);
    catch
        values(k) = NaN;
    end
end

% Update akamat
for k = 1:nPending
    i = pending_i(k);
    j = pending_j(k);

    ii = j - pmin + 1;
    akamat(ii, i+1) = values(k);
end
total_time = toc;

% Verbose output and timing summary
if verbflag
    fprintf('Completed in %.f seconds.\n', total_time);
end

% File writing (if required)
if ~isempty(iw)
    % Write the entire matrix to the file
    dlmwrite(nomfich, akamat, 'delimiter', ' ', '-append'); %#ok<DLMWT>
end

% Prepare outputs
varargout = {akamat, pmax, qmax, N};
end