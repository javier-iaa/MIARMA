function plotint(s1, s2, interp, yfor, yback)
% function plotint(s1, s2, interp) plots the interpolated segment between
% s1 and s2 for test purposes.
% Optional arguments:
% varargin{1} is the forward extrapolation
% varargin{2} is the backward extrapolation
%
% Version: 0.1
% Changes from the last version:
%
%  Author(s): Javier Pascual-Granado
%  Date: 19/03/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure;
    ls1 = length(s1);
    ls2 = length(s2);
    np = length(interp);

    % Plot data segments
    plot(1:ls1, s1, 'b.')
    hold
    plot((ls1+np+1):(ls1+np+ls2), s2, 'm.')

    % Plot the interpolation
    plot((1+ls1):(ls1+np), interp, 'k.')

    % Plot extrapolations
    if nargin > 3
        plot(ls1:(ls1+np), yfor, 'r.')
        plot((1+ls1):(1+ls1+np), yback, 'g.')
    end

end