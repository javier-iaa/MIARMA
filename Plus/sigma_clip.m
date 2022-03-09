function [ seg_out, stat_out ] = sigma_clip(seg_in, varargin)
% function seg_out = sigma_clip(seg_in, varargin) performs a sigma clipping
% of the input data segment seg_in with a sigma factor fac (default is 2.5).
% seg_in is initially detrended with a polynomial of 2nd degree.
% The output data segment seg_out substitute values where seg_in 
% exceed the limit with the mean of detrended segment.
% 
% Optional inputs are: 'fac' for sigma factor, 'verb' for verbosity and 'empty' to
% fill clipped data with zeros.
% stat_out is by default a vector of zeros but if the 'empty' flag is used clipped 
% values  are set to ones
%
% Version: 0.4.1
% Changes from the last version:
% - Info is updated and fixed.
%
%  Author(s): Javier Pascual-Granado
%  Date: 12/01/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Verbosity flag
    verb_flag = find( strcmp(varargin,'verb'), 1 );
    if ~isempty(verb_flag)
        verb = true;
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
    
    % leave empty flag
    empty_flag = find( strcmp(varargin, 'empty'), 1 );
    
    l1 = numel(seg_in);
    
     % status vector
    stat_out = ones(1, l1);
    
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
    if ~isempty(empty_flag)
        seg_out( iclip ) = 0;
        stat_out( iclip ) = 1;
    else
        y = seg_det( inclip );
        sclip = interp1( inclip, y, iclip, 'pchip');
        seg_out( iclip ) = sclip + f( iclip );
    end
    
    % Reshaping into cols.
    seg_out = reshape(seg_out, l1, 1);
    
end
