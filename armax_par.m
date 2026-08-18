function model1 = armax_par(z, o)
% function model1 = armax_par(z, o) model the data segment z by using an
% arma loss function which is minimised through fmincon.
% 
% The algorithm used for minimisation, SQP, can be parallelised so this
% might be more effective than armax function when many cores are
% available. Otherwise, it is more optimal to use armax.
%
% Inputs:       z - data segment to be modelled
%               o - ARMA (p,q) orders
%
% Outputs:      model1 - ARMA model
%
% Version: 0.1
%
%  Author(s): Javier Pascual-Granado
%  Date: 17/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs: supply your iddata object z, and model orders o
na = o(1); nc = o(2);   % ARMA orders

% Seed (optional)
% try
%     ord_armax = [na 0 nc 0];
%     seed = armax(z, ord_armax);
%     p0 = [seed.A(2:end), seed.C(2:end)]; % drop leading 1
% catch
%     p0 = 0.1*ones(1, na+nc);
% end

p0 = 0.1*ones(1, na+nc);

% Objective for ARMA (SISO)
function J = arma_loss(p, z, na, nc)
    if any(~isfinite(p))
        J = 1e20;
        return
    end
    A = [1, p(1:na)];
    C = [1, p(na+1:na+nc)];
    sys = idpoly(A, [], C);             % no B polynomial
    if ~isstable(sys)
        J = 1e20;
        return
    end
    try
        yhat = predict(sys, z, 1);          % one-step-ahead prediction
    catch
        J = 1e20;
        return
    end
    y = z;
    err = y - yhat;
    J = sqrt(mean(err.^2)) / std(y);    % normalized RMS (like loss)
    J = real(J);
end

% Optimizer options (parallel finite-difference)
opts = optimoptions('fmincon', ...
     'Display','none', ...
     'Algorithm','sqp', ...
     'UseParallel', true, ...
     'MaxIterations', 200);

% Optional: suppress warnings temporarily
warnState = warning('off','all');

lb = -10*ones(1, na+nc);
ub =  10*ones(1, na+nc);

if isempty(gcp('nocreate'))
    % start pool for parallel computing without printing
    evalc('parpool'); 
end

obj = @(p) arma_loss(p, z, na, nc);
[popt, Jmin] = fmincon(obj, p0, [], [], [], [], lb, ub, [], opts);

% Reconstruct final ARMA model
Aopt = [1, popt(1:na)];
Copt = [1, popt(na+1:na+nc)];
model1 = idpoly(Aopt, [], Copt);

% Print the loss function
% fprintf('Optimized loss = %.6g\n', Jmin);

% Restore warning state
warning(warnState);

end