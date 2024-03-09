function interp = polintre(seg1, seg2, np, ord)
% Interpolates between segments seg1 and seg2, np datapoints using a
% polynomial fit with order ord. Estimate statistical properties of the
% residuals and add a stochastic component mimicking the noise.
%
% Note: polintre is prepared to make forward and backward extrapolation too
% by design.
%
% Version: 0.1
% Changes from the last version:
% - Optimization and minor fixes.
% - Note the factor 0.1 in line 46 to limit noise dispersion.
%
% Author(s): Javier Pascual-Granado
% Date: 04/03/2021
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare data
[fil,~] = size(seg1);
if fil==1
    seg1 = seg1';
end
[fil,~] = size(seg2);
if fil==1
    seg2 = seg2';
end
data = [seg1; seg2];
l1 = length(seg1);
l2 = length(seg2);
L = l1 + np + l2;
t1 = 1:l1;
t_int = (l1+1):(l1+np);
t2 = (l1+np+1):(l1+np+l2);
t_tot = 1:L;
t_data = [t1 t2];

% Fit model
pp_data = polyfit(t_data', data, ord);
pval_tot = polyval(pp_data, t_tot);
pval_int = pval_tot(t_int);

% Get residuals
res = data - pval_tot(t_data)' ;
rstd = std(res);
r = 0.1*rstd*randn(size(t_int));

% Interpolation
interp = pval_int + r;
