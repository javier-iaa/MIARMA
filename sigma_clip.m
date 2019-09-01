function seg_out = sigma_clip(seg_in,fac)
% function seg_out = sigma_clip(seg_in,fac) performs a sigma clipping
% of the input data segment seg_in with a sigma factor fac.
% seg_in is initially detrended with a polynomial of 2nd degree.
% The output data segment seg_out substitute values where seg_in 
% exceed the limit with the mean of detrended segment.

    l1 = numel(seg_in);
    
    % Detrending
    x = 1:l1;
    
    % Reshaping into rows
    seg_in = reshape(seg_in,1,l1);
    
    % Detrending with a 2-polynomial
    p = polyfit(x,seg_in,2);
    f = polyval(p,x);
    seg_out = seg_in - f;
    
    % Sigma clipping
    med = mean(seg_out);
    sig = std(seg_out);
    ind = abs(seg_out-med) > fac*sig;
    seg_out(ind) = med;
    seg_out(~ind) = seg_out(~ind) + f(~ind);
    
    % Reshaping into cols.
    seg_out = reshape(seg_out,l1,1);
end
