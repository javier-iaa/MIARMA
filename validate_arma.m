function [isval, sta, interp] = validate_arma( data, ord, facint, check, memtot, interp )
% Function isval = validate_arma( data, ord, facint, check ) validate the
% ARMA model of data <data> with order <ord> using the test <check>.
%
% In order to perform the diagnostics a single gap is introduced in the middle 
% of the time series and armaint.m is used to interpolate. Then, several 
% statistics are used to compare between original and interpolated data.
% 
% The gap introduced will have a size
%
%   gapsize = L / (2*facint + 1)
%
% where L is the length of data and <facint> is an input from MIARMA.
%
% In this way, each data segment bracketing the gap will have the same 
% length,  which is exactly a factor facint times the gap size.
%
% Several diagnostics can be used to test the interpolation and validate 
% ord as the optimal order. Parameter <check> can be any of:
% 'fvar' for fraction of variability
% 'pc' for Percentage of Consistency (default)
% 'wn' for whiteness test
% 'mse' for Mean Squared Error
%
% The output isval is true if the test is passed and sta is the value of 
% the corresponding statistic used.
% 
% By Javier Pascual-Granado
% <a href="matlab:web http://www.iaa.es;">IAA-CSIC, Spain</a>
%
% Version: 0.2.3
%
% Changes: 
% - Added limit in segment length to avoid memory overflow when the number
% of free parameters of the model is very high.
%
% Call: armaint.m
%
% Date: 12/09/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = length(data);

if ~exist('check','var')
    check = 'pc';
end

% Convert to column-wise vector
data = reshape(data, L, 1);

% Statistical normalization
sn = ( data - mean(data) )/std(data);

% Data segments
gapsize = fix( L / (2*facint + 1) );
segsize = facint*gapsize;

po = ord(1);

% This is approx. the segment limit to avoid an overflow
lim_segsize = floor( memtot/(8*(po^3 + 4*po^2 + po)/1024^3) - 10 );

if segsize > lim_segsize
    segsize = lim_segsize;
end

indf1 = segsize;
seg1 = sn(1:indf1);
orig = sn( (indf1+1):(indf1+gapsize) );
indi2 = indf1 + gapsize + 1;
indf2 = L;
seg2 = sn(indi2:indf2);

if length(seg2) > lim_segsize
    seg2 = seg2(1:lim_segsize);
end

% Interpolation
if ~exist('interp', 'var')
    [interp, ~] = armaint(seg1, seg2, ord, gapsize);
end

% Mean of the data segment
m = mean(orig);

if any( isnan(interp) )
    isval = false;
    sta = NaN;
    return
end

% Prediction error
res = interp - orig;

%% 1st check: compare variability of model and measurements
% Percentage of output variability that is explained by the model. This is 
% similar to the fit parameter calculated by Matlab compare function
if strcmp( check, 'fvar')
    
    lim_fvar = 50;

    sta = 100*( 1 - norm(res)/norm(orig-m) ) ;

    if sta<lim_fvar
        isval = false;
    else
        isval = true;
    end
end
%% 2nd check: Consistency check (Ding et al. 2000)
if strcmp( check, 'pc')
     % Auto-correlation of original data
    [c_orig,~] = xcorr(orig, floor(gapsize/2), 'coeff');

    % Auto-correlation of interpolated data
    [c_interp,~] = xcorr(interp, floor(gapsize/2), 'coeff');

    % Percentage of Consistency (>80% or 85% is ok)
    sta = 100*( 1-sum( (c_interp - c_orig).^2 ) / sum(c_orig.^2) );

     if sta<80 || isnan(sta)
         isval = false;
     else
         isval = true;
     end
 end
%% 3rd check: whiteness
% Not implemented yet
if strcmp( check, 'wn'  )
    m1 = armax(seg1, ord);
    data1 = iddata( seg1, []);
    resid(m1, data1);
end
 
%% 4th check: Mean Squared Error (MSE)
% Not working fine

% lim_chi2 = 0.2;
% 
% if strcmp( check, 'mse')
%     samp = 60;
%     % Use abs if orig has zero-mean
%     squn = sqrt( abs(orig)/samp ) .^2;
%     sres = (orig - interp').^2;
%     sta = sum( sres./ squn ) / gapsize;
% 
%     if abs(sta-1)>lim_chi2 || isnan(sta)
%         isval = false;
%     else
%         isval = true;
%     end
% end

if strcmp( check, 'mse')
    res = abs(interp - orig);
    [h, sta] = chi2gof( res );
    
    if h || isnan(sta)
        % The null hypothesis cannot be rejected (i.e. correlated residuals) 
        % or there are NaN in the residuals
        isval = false;
    else
        isval = true;
    end
end

end
