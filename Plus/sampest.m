function samp = sampest(timein)
% function samp = sampest(timein) estimates the sampling of the time series
% with times in timein.
%
%  Version: 0.3.1
%  Changes: histogram of the sampling intervals. This will be deprecated in
%  future versions using the Bayesian Blocks algorithm to calculate the
%  histogram since the classical histogram gives a poor resolution and it
%  is necessary a huge number of bins.
%
% Changes: changed mean into median and simplified. Samp0 and err0 fixed.
% Number of bins adjusted. When histogram approach fails uses median
% instead.
%
%  Author: Javier Pascual-Granado
%  $Date: 20/10/2023 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


L = length(timein);

% Time interval between measures
dt = diff( timein);

%% Estimate sampling

sampm = median(dt);

N = L-1; % number of bins

[n, dtm] = hist(dt, N);
[~, I] = max(n);
samp0 = 2*dtm(I);

err0 = samp0/2;

% final estimate is obtained
samp = median( dt( abs(dt-samp0)<err0 ) );

if isnan(samp) | sampm<samp
    samp = sampm;
end
