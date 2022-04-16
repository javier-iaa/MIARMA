function [av_SXX, av_SXXh, Ndata] = fastCGSA(data, h, l, sub) 
%[av_SXX, av_SXXh, Ndata] = fastCGSA(data , h, l, sub) 
%
% Input:
%              data is the original time series
%
%              h is the scale for coarse graining
%               default is 2.
%
%              l is the fraction of data to use for subsets
%               e.g. l=0.9  Yamamoto; l=0.25 Jung Shao
%
%              sub is the number of subsets.
%               default is 10.
%
% Output:
%               av_SXX is the auto-power spectrum averaged through the sub 
%               subsets.
% 
%               av_SXXh is the cross-power spectrum between the time series
%               X and the coarse grained version Xh averaged through the
%               sub subsets.
%
%               Ndata is the number of data points used from data
%
% Version: 0.1
%
% Code based on Yamamoto and Hughson 1991 and 1993.
%
% Author: Javier Pascual-Granado
% Last update: 08-sep-2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Default values
if nargin==1,
    h = 2;
    l = 0.9;
    sub = 10;
end

L = length(data);

% The length is truncated to the nth in order to use the h scale
Ntot = floor(L/h);
X = data(1:Ntot);

% convert columns into rows
X = reshape(X,1,Ntot);

% Calculate subsets length
Ndata = 2^floor( log2( l*Ntot ) );

% Calculate lag between subsets
lag = floor( (Ntot - Ndata) / (sub - 1) );

% Coarse grained time series
X_h = data(1:h:end);

% Define the cosine-taper window
w = tukeywin(Ndata, 0.5);

% Initialize the DFT matrices
FS = zeros(sub, Ndata);
FS_h = zeros(sub, Ndata);

% Taper and FFT segments of data
for i=0:sub-1
    % Segment indexes
    indi = 1+i*lag;
    indf = Ndata+i*lag;
    ind = indi:indf;
    
    % Tapering of each segment
    S = X(ind) .* w';
    Sh = X_h(ind) .* w';
 
    % DFT of each segment
    FS(i+1,:) = fft(S)/Ndata;
    FS_h(i+1,:) = fft(Sh)/Ndata;
end

% Averaged auto-power spectrum
av_SXX = (1/sub)*sum( abs( FS(3:end, :) ).^2 , 1);

% Cross-power spectra
csp_h = FS.*conj(FS_h);

% Averaged through the subsets
av_SXXh = (1/sub)*sum( abs( csp_h ) , 1);

end
