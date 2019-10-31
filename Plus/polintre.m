function interp = polintre(seg1, seg2, np, ord)
% Interpolates between segments seg1 and seg2, np datapoints using a
% polynomial fit with order ord. Estimate statistical properties of the
% residuals and add a stochastic component mimicking the noise.

% Weights for the interpolation
% wp = 1/(np+1);
% w = wp:wp:(1-wp)

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
pval_seg1 = pval_tot(t1);
pval_seg2 = pval_tot(t2);
pval_int = pval_tot(t_int);

% Get residuals
res1 = seg1 - pval_seg1';
res2 = seg2 - pval_seg2';
res = [res1; res2];
rstd = std(res);
r = rstd*randn(size(t_int));

% Interpolation
interp = pval_int + r;