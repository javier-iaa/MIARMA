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
% 'rep' followed by an integer sets the limit in the repetition of orders, 
% otherwise repetition is not used as a stop criterion.
%
% By Javier Pascual-Granado
% <a href="matlab:web http://www.iaa.es;">IAA-CSIC, Spain</a>
%
% Version: 0.2.2 R2022
%
% Changes:
% - Reuses interp from previous validate_arma call.
%
% Calls:
% validate_arma 0.2.2
%
% Date: 08/07/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

verbflag = true;
temp = true;
ML = length(seg);

% Check inputs
optargin = size(varargin, 2);

if optargin==0
    temp = false;
elseif optargin==1
    if ~strcmp(varargin{1},'w')
        error('Wrong inputs');
    end
else
    iw = find(strcmp(varargin, 'w'), 1);
    if ~isempty(iw)
        if ischar(varargin{iw+1})
            akaname = varargin{iw+1};
        end
    end
    
    iv = find(strcmp(varargin, 'verbose'), 1);
    if ~isempty(iv)
        if islogical(varargin{iv+1})
            verbflag = varargin{iv+1};
        else
            error('Wrong inputs');
        end
    end
    
    ir = find(strcmp(varargin, 'rep'), 1);
    if ~isempty(ir)
        rep_lim = varargin{ir+1};
    else
        rep_lim = 10;    % default value for repetition 
    end
    
    im = find(strcmp(varargin, 'mseg'), 1);
    if ~isempty(im)
        % mseg is the max. length of the segment use to calculate the order
        mseg = varargin{im+1};
        
        % If the statistical tests are not passed Akaike matrix will be recalculated 
        % using seg0 instead of seg
        seg0 = seg;
        if ML>mseg
            seg = seg(1:mseg);
        end

        if verbflag
            fprintf('\n%d points will be used for finding optimal order.\n', length(seg));
        end
    end
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
lseg = length(seg);

while lseg<=ML

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

        % Stop condition: if the same order is found in more than rep_lim iterations
        if exist('rep_lim', 'var')
            if p==p0 && q==q0
                cont = cont + 1;
                if cont==rep_lim
                    break
                end
            else
                cont = 0;
            end
            p0 = p;
            q0 = q;
        end

        % Validate model with optimal order
        [isval_mse, sta_mse, interp] = validate_arma( seg, ord, facint, 'mse' );
        [isval_pc, sta_pc] = validate_arma( seg, ord, facint, 'pc', interp );
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
            break
        end

    end
    
    if isval
        return
    end
    
    if lseg<ML
        seg = seg0;
        lseg = ML;
    else
        return
    end
    
end
