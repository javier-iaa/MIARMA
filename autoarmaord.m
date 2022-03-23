function aka = autoarmaord( seg, varargin)
% function aka = autoarmaord( seg, varargin)
% Call armaord.m iteratively until a model with an 
% optimal order is validated.
% 
% Possible calls include:
% aka = autoarmaord( seg );
% No Akaike file is saved.
% Verbosity is on.
% 
% aka = autoarmaord( seg, 'w' );
% Akaike matrix is saved in file with default name.
% Verbosity is on.
% 
% aka = autoarmaord( seg, 'w', akaname);
% Akaike matrix is saved in file with name <akaname>.
% Verbosity is on.
% 
% aka = autoarmaord( seg, 'verbose', false );
% No Akaike file is saved. 
% Verbosity is set off.
% 
% aka = autoarmaord( seg, 'verbose', false, 'w');
% Akaike matrix is saved in file with default name.
% Verbosity is set off.
% 
% aka = autoarmaord( seg, 'verbose', false, 'w', akaname);
% Akaike matrix is saved in file with name <akaname>.
% Verbosity is set off.
% 
% By Javier Pascual-Granado
% <a href="matlab:web http://www.iaa.es;">IAA-CSIC, Spain</a>
%
% Version: 0.1
%
% Changes:
% - Introduced stop condition when the optimal order repeats several times
% - Minor fixes.
%
% Calls:
%  - armaord.m
%  - validate_arma.m
%
% Date: 22/03/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check inputs

optargin = size(varargin, 2);
verbflag = true;
temp = true;

if optargin==0
    temp = false;
elseif optargin==1
    if ~strcmp(varargin{1},'w')
        error('Wrong inputs');
    end
elseif optargin==2
    if strcmp(varargin{1},'w')
        if strcmp(varargin{2},'verbose')
            error('Wrong inputs');
        else
            akaname = varargin{2};
        end
    elseif strcmp(varargin{2},'w')
            error('Wrong inputs');
    elseif strcmp(varargin{1}, 'verbose')
        verbflag = varargin{2};
    end
elseif optargin==3
    if ~strcmp(varargin{end},'w')
        error('Wrong inputs');
    elseif strcmp(varargin{1},'verbose')
        verbflag = varargin{2};
    else
        error('Wrong inputs');
    end
elseif optargin==4
    if ~strcmp(varargin{3},'w')
        error('Wrong inputs');
    elseif strcmp(varargin{1}, 'verbose')
        akaname = varargin{4};
        verbflag = varargin{2};
    else
        error('Wrong inputs');
    end
elseif optargin>4
    error('Wrong number of inputs');
end

% Initial conditions
isval = 0;
pmin = 2; % pmin, facint will be inputs in next update
facint = 1.5;
dpmax = 16; % dpmax, dqmax might be optional 
dqmax = 6;   % inputs too in next update
pmax = dpmax;
qmax = dqmax;
pmax_lim = dpmax*10; % pmax_lim, qmax_lim might be 
qmax_lim = dqmax*10; % optional inputs too in next update
p0 = 0;
q0 = 0;
cont = 0;

while ~isval
    
    % Calculate Akaike matrix
    if exist( 'akaname', 'var' )
        if verbflag
            aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
                'qmax', qmax, 'w', akaname);
        else
            aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
            'qmax', qmax, 'verbose', false, 'w', akaname);
        end
        
    elseif temp
        if verbflag
            aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
                'qmax', qmax, 'w' );
        else
            aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
            'qmax', qmax, 'verbose', false, 'w' );
        end
    else
        if verbflag
            aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
                'qmax', qmax);
        else
            aka = armaord( seg, 'pmin', pmin, 'pmax', pmax, ...
            'qmax', qmax, 'verbose', false);
        end
    end
    
    % Find optimal order
    [cp, cq] = find( aka == min( min( aka ) ) );
    q = cq - 1;
    p = cp + pmin - 1;
    ord = [p q];
    fprintf('\nOptimal order is (%d, %d)\n', p, q);
    
    % Stop condition: if the same order is found in more than 2 iterations
    if p==p0 && q==q0
        cont = cont + 1;
        if cont==2
            return
        end
    else
        cont = 0;
    end
    p0 = p;
    q0 = q;

    % Validate model with optimal order
    [isval_mse, sta_mse] = validate_arma( seg, ord, facint, 'mse' );
    [isval_pc, sta_pc] = validate_arma( seg, ord, facint, 'pc' );
    if isval_mse
        str_mse = 'cannot reject';
    else
        str_mse = 'reject';
    end
    if isval_pc
        str_pc = 'cannot reject';
    else
        str_pc = 'reject';
    end
    fprintf('\nMSE test %s null hypothesis with h=%4.1e\n', str_mse, sta_mse);
    fprintf('PC test %s null hypothesis with h=%4.1e\n', str_pc, sta_pc);
    isval = isval_mse && isval_pc;
    
    pmax = pmax + dpmax;
    qmax = qmax + dqmax;
    
    if pmax>pmax_lim || qmax>qmax_lim
        return
    end
    
end
