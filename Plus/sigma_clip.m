function seg_out = sigma_clip(seg_in, varargin)
% function seg_out = sigma_clip(seg_in, varargin) performs a sigma clipping
% of the input data segment seg_in with a sigma factor fac.
% seg_in is initially detrended with a polynomial of 2nd degree.
% The output data segment seg_out substitute values where seg_in 
% exceed the limit with the mean of detrended segment.
% 
% Optional inputs are: 'fac', 'verb'
%
% Version: 0.3
% Changes from the last version:
%  - Using pchip interpolation instead of the mean
%  - few improvements and bug fixes.
%
%  Author(s): Javier Pascual-Granado
%  Date: 29/03/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Verbosity flag
    verb_flag = find( strcmp(varargin,'verb') );
    if ~isempty(verb_flag)
        verb = varargin{verb_flag + 1};
    else
        verb = false; % Default is off
    end

    % factor flag
    fac_flag = find( strcmp(varargin, 'fac') );
    if ~isempty(fac_flag)
        fac = varargin{fac_flag + 1};
    else
        fac = 2.5; % Default
    end
 
    l1 = numel(seg_in);
    
    % Detrending
    x = 1:l1;
    
    % Reshaping into rows
    seg_in = reshape(seg_in,1,l1);
    seg_out = seg_in;
    
    % Detrending with a 2-polynomial
    p = polyfit(x,seg_in,2);
    f = polyval(p,x);
    seg_det = seg_in - f;
    
    % Sigma clipping
    med = mean(seg_det);
    sig = std(seg_det);
    ind = abs(seg_det-med) > fac*sig;
    iclip = find( ind );
    inclip = find( ~ind );
    
    if verb
        nc = length( iclip );
        fprintf('%d clipped points\n', nc);
    end
    
    % clipped datapoints
    y = seg_det( inclip );
    sclip = interp1( inclip, y, iclip, 'pchip');
    seg_out( iclip ) = sclip + f( iclip );
    
    % Reshaping into cols.
    seg_out = reshape(seg_out, l1, 1);
    
end
