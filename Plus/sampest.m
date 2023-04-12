function samp = sampest(timein)
% function samp = sampest(timein) estimates the sampling of the time series
% with times in timein.
%
%  Version: 0.3
%  Changes: histogram of the sampling intervals. This will be deprecated in
%  future versions using the Bayesian Blocks algorithm to calculate the
%  histogram since the classical histogram gives a poor resolution and it
%  is necessary a huge number of bins.
%
% Changes: changed mean into median and simplified.
%
%  Author: Javier Pascual-Granado
%  $Date: 12/04/2023 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


L = length(timein);

% Time interval between measures
dt = diff( timein);

%% Estimate sampling
[n, dtm] = hist(dt, L-1);
[~, I] = max(n);
samp0 = dtm(I);

err0 = samp0/2; 

% final estimate is obtained
samp = median( dt( abs(dt-samp0)<err0 ) );                                          
