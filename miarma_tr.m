function miarma_tr( fname, transit_str, varargin)
% This is a wrapper to run MIARMA on a file removing transits before the 
% gap-filling algorithm is applied.
% 
% Apart from the filename string <fname> a structure <transit_str> must be 
% provided with the following fields:
%    transit_str.porb       -> orbital period in days
%    transit_str.epoch     -> mid-transit in RJD
%    transit_str.duration -> transit duration in days
%
% Other inputs that might be provided preceded by the corresponding tagname 
% string are:
%   'nhead' is the number of header lines in the ASCII file (default 1)
%   'delim' is the string delimiter to use (default is ' ')
%   'timecol' is the column number for time (default is 1)
%   'magcol' is the column number for magnitude (default is 2)
%   'statcol' is the column number for the status of datapoints, i.e. the flag 
%     to decide when to interpolate or not. Status is modified after transit removal
% 
% Example for Kepler data
% 
% fname = 'kplr007199397-2011116030358_slc.dat';
% transit_str.porb = 105.881767; 
% transit_str.epoch = 2454989.979350 - 2400000.0; 
% transit_str.duration = 18.0440;
% 
% miarma_tr(fname, transit_str, 'nhead', 8, 'delim', ' ', 'timecol', 1, 'magcol', 4);
% 
% 
% By Javier Pascual-Granado
% <a href="matlab:web http://www.iaa.csic.es;">IAA-CSIC, Spain</a>
%
% Version: 0.1.2 - R2022
% Changes:
% - Minor fixes.
%
% Date: 03/06/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Import data
porb = transit_str.porb;
epoch = transit_str.epoch;
duration = transit_str.duration;

% Value for delim parameter
idl = find(strcmp(varargin,'delim'), 1);
if isempty(idl)
    delim = ' '; % Default
else
    delim = varargin{idl+1};
end

% Value for nhead parameter
inh = find(strcmp(varargin,'nhead'), 1);
if isempty(inh)
    nhead = 1; % Default
else
    nhead = varargin{inh+1};
end

% Value for timecol parameter
itc = find(strcmp(varargin,'timecol'), 1);
if isempty(itc)
    timecol = 1; % Default
else
    timecol = varargin{itc+1};
end

% Value for magcol parameter
imc = find(strcmp(varargin,'magcol'), 1);
if isempty(imc)
    magcol = 2; % Default
else
    magcol = varargin{imc+1};
end

% Import data
data = importdata(fname, delim, nhead);
data = data.data;

% Time
t = data(:, timecol);
L = length(t);

% Magnitude
s = data(:, magcol);

% Value for statcol parameter
isc = find(strcmp(varargin,'statcol'), 1);
if isempty(isc)
    stat = zeros(1,L); % Default
else
    statcol = varargin{isc+1};
    stat = data(:, statcol);
end

% Remove inf is there is any
indinf = find( isinf(s) );
s(indinf) = [];
t(indinf) = [];
stat(indinf) = [];
%linf = length( indinf );

% Convert to ppm and detrend
% p = polyfit(t, s,1);
% y = polyval(p, t);
% sd = 1e6*(s./y - 1);

%% Remove transits

% Transit removal
trans_ind_ini = find( mod(t-epoch, porb)<0.001,1 ); % first mid-transit index

if isempty(trans_ind_ini)
    fprintf('No transit found\n');
    return
end
transitt = t(trans_ind_ini):porb:t(end); % transit times
N = length(transitt);

sc = s;
tc = t;
statc = stat;

for i=1:N
    ttini = transitt(i) - duration/2;
    ttfin = transitt(i) + duration/2;
    transit_ind = find(tc > ttini & tc<ttfin);
    tc(transit_ind) = [];
    sc(transit_ind) = [];
    statc(transit_ind) = [];
end

%lost_transit = L - length(tc);

%% MIARMA

% First-step: sampling regularization (no longer necessary since it is implemented inside MIARMA)
% [treg, sreg, statreg] = regsamp(tc, sc, statc);

% Inputs for MIARMA
% strdata.time = treg;
strdata.time = tc;
% strdata.data = sreg;
strdata.data = sc;
% strdata.stat = statreg;
strdata.stat = statc;

% strdata.params.temp = 1;
% strdata.params.pmax = 20;
% strdata.params.facmin = 10;
strdata.params.mseg = 3000;

% Gap-filling
strout = MIARMA( strdata );
sint = strout.datout;
tint = strout.timeout;

%% Write the output
Lg = length(sint);
fnameo = [ fname(1:end-4) '_miarma.dat' ];
fich = fopen(fnameo, 'w');
for i=1:Lg
    %fprintf(fich,'%16.13f %16.13f\n', treg(i), sint(i));
    fprintf(fich,'%16.13f %16.13f\n', tint(i), sint(i));
end
fclose(fich);
