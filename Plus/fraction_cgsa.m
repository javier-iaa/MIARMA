function frac_per = fraction_cgsa(X)
% This function calculates the fraction of fractal component present in 
%  the time series X using CGSA algorithm

    [x, xh, ~] = fastCGSA(X);
     
    % Fractal percentage
    frac_per = 100*sum(xh)/sum(x);

end
