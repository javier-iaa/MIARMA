function plotint(s1, s2, interp, yfor, yback)
% function plotint(s1, s2, interp, yfor, yback) plots the interpolated segment between
% s1 and s2 for test purposes.
% Optional arguments:
% yfor and yback are forward and backward extrapolations respect.
%
% s1 in blue
% interp in black
% s2 in magenta
% yfor in red
% yback in green
%
% Version: 0.2
% Changes from the last version:
% - Bug fixes
% - Marks increased.
%
%  Author(s): Javier Pascual-Granado
%  Date: 20/05/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure;
    ls1 = length(s1);
    ls2 = length(s2);
    np = length(interp);
    
    if isnan(interp)
        fprintf('interp is nan');
        return
    end

    % Plot data segments
    plot(1:ls1, s1, 'b.', 'MarkerSize', 18)
    hold
    plot((ls1+np+1):(ls1+np+ls2), s2, 'm.', 'MarkerSize', 18)

    % Plot the interpolation
    plot((1+ls1):(ls1+np), interp, 'k.', 'MarkerSize', 18)

    % Plot extrapolations
    if nargin > 3
        plot(ls1:(ls1+np), yfor, 'r.', 'MarkerSize', 18)
        plot((1+ls1):(1+ls1+np), yback, 'g.', 'MarkerSize', 18)
    end

end
