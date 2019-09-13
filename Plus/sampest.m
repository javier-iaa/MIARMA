function samp = sampest(timein)
% function samp = sampest(timein) estimates the sampling of the time series
% with times in timein.
%
%  Version: 0.1
%  Changes: histogram of the sampling intervals. This will be deprecated in
%  future versions using the Bayesian Blocks algorithm to calculate the
%  histogram since the classical histogram gives a poor resolution and it
%  is necessary a huge number of bins.
%
%  Author: Javier Pascual-Granado
%  $Date: 02/03/2019 $
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


L = length(timein);

% Time interval between measures
dt = timein(2:L)-timein(1:L-1);

%% Estimate sampling
% a 1-percent deviation is admitted
[n, dtm] = hist(dt,10000);
[~, I] = max(n);
samp0 = dtm(I);
% samp0 = timein(2)-timein(1);
err0 = samp0/2; 
% ijump = find(abs(dt-samp0)>err0,1);

% the estimate is refined
% samp1 = mean(dt(1:(ijump-1)));
% err1 = samp1/100;

% final estimate is obtained
samp = mean( dt( abs(dt-samp0)<err0 ) );
% err = samp/100;                                                 
