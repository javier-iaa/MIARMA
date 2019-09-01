function varargout = regsamp(varargin)
% function varargout = regsamp(varargin) creates a regularly sampled
% version of a time series introducing time measures inside gaps with a
% given regular sampling.
% An approximately regular sampling is assumed in the segments without
% gaps.
% Inputs:   an ASCII filename of a file containing the input data array, the
%           corresponding time array, and the status array
%
%           For test purposes this call is also available:
%           regsamp(timein,datain,flagin)
%
% Outputs:  If no output variable is given an ASCII file is created 
%           containing the output data and the corresponding time and
%           status, otherwise:
%           [timeout,dataout,flagout] = regsamp(timein,datain,flagin)
%
% Version: 0.6
% Changes: Time values are extrapolated to get evenly spaced time series
%
%  Author: Javier Pascual-Granado
%  $Date: 14/03/2019$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Input data
if ischar(varargin{1}),
    filename = varargin{1};
    strdata = importdata(filename);
    timein = strdata.data(:,1);
    datin = strdata.data(:,2);
    flagin = strdata.data(:,3);
else
    timein = varargin{1};
    datin = varargin{2};
    flagin = varargin{3};
end

L = length(timein);

% Change columns into rows
timein = reshape(timein,1,L);
datin = reshape(datin,1,L);
flagin = reshape(flagin,1,L);

% Time interval between measures
dt = timein(2:L)-timein(1:L-1);

%% Estimate sampling
samp = sampest(timein); % estimate of the mean global sampling
err = samp/1000;
ind = find((dt-samp)>err);

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

timein = timein - timein(1);

timef = samp*round(timein(end)/samp);
timeout = 0:samp:timef;
L1 = length(timeout);
datout = zeros(1,L1);
flagout = ones(1,L1);
i1 = zeros(1,L); % indexes of the nearest timeout to timein values

for i=1:L,
    [~, i1(i)] = min(abs(timeout-timein(i))); 
end

new_do = interp1(timein, datin, timeout(i1),'linear','extrap');
datout(i1) = new_do;
flagout(i1) = flagin;

varargout{1} = timeout;
varargout{2} = datout;
varargout{3} = flagout;

end