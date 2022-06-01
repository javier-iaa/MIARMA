function varargout = regsamp(varargin)
% function varargout = regsamp(varargin) creates a regularly sampled
% version of a time series introducing time measures inside gaps with a
% given regular sampling.
% An approximately regular sampling is assumed in the segments without
% gaps.
% In any case, time values are extrapolated to get evenly spaced time series
%
% Inputs:   an ASCII filename of a file containing the input data array, the
%           corresponding time array, and the status array (optional)
%
%           For test purposes this call is also available:
%           regsamp(timein,datain,flagin)
%
% Outputs:  If no output variable is given an ASCII file is created 
%           containing the output data and the corresponding time and
%           status, otherwise:
%           [timeout,dataout,flagout] = regsamp(timein,datain,flagin)
%
% Version: 0.8.3
% Changes: 
% - Using median values instead of zeros.
%
% Requires ismemberf.m from Bruno Luong
%
%  Author: Javier Pascual-Granado
%  $Date: 14/04/2022$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input data
if ischar(varargin{1}),
    filename = varargin{1};
    strdata = importdata(filename);
    timein = strdata.data(:,1);
    datin = strdata.data(:,2);
    if size(strdata.data,2)==3
        flagin = strdata.data(:,3);
    else
        flagin = zeros( size(timein) );
    end
else
    timein = varargin{1};
    datin = varargin{2};
    if nargin==3
        flagin = varargin{3};
    else
        flagin = zeros( size(timein) );
    end
end

L = length(timein);

% Change columns into rows
timein = reshape(timein,1,L);
datin = reshape(datin,1,L);
flagin = reshape(flagin,1,L);

% Time interval between measures
dt = diff( timein );

%% Estimate sampling
samp = sampest(timein); % estimate of the mean global sampling
err = samp/1000;
ind = find((dt-samp)>err,1);

if isempty(ind),
    if ischar(varargin{1}),
        fout = cat(2,filename(1:end-4),'.rgs');
        system(cat(2,'cp ',filename,' ',fout));
    else
        varargout{1} = timein;
        varargout{2} = datin;
        varargout{3} = flagin;
    end
    return;
end

 t0 = timein(1);
 timein = timein - t0;

timef = samp*round( timein(end)/samp );
timeout = 0:samp:timef;
L1 = length(timeout);
med = median( datin );
% datout = zeros(1,L1);
datout = med*ones(1,L1);
flagout = ones(1,L1);

ti = samp*nearest( timein/samp );    
errti = 10*eps( ti(end) );

% indexes of the nearest timeout to timein values
[~, i1] = ismemberf(ti, timeout, 'tol', errti);

new_do = interp1( timein, datin, timeout(i1), 'pchip', 'extrap' );
datout(i1) = new_do;

flagout(i1) = flagin;
timeout = timeout + t0;

varargout{1} = timeout;
varargout{2} = datout;
varargout{3} = flagout;

end
