function model1 = armax_par(z, o)
% Inputs: supply your iddata object and model orders
%z = seg1n;                     % iddata
na = o(1); nc = o(2);   % ARMA orders
ord_armax = [na 0 nc 0];   % for seed estimation

% Seed (optional)
try
    seed = armax(z, ord_armax);
    p0 = [seed.A(2:end), seed.C(2:end)]; % drop leading 1
catch
    p0 = 0.1*ones(1, na+nc);
end

% Objective for ARMA (SISO)
function J = arma_loss(p, z, na, nc)
    A = [1, p(1:na)];
    C = [1, p(na+1:na+nc)];
    sys = idpoly(A, [], C);             % no B polynomial
    yhat = predict(sys, z, 1);          % one-step-ahead prediction
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

opts = algoprop();

% Optional: suppress warnings temporarily
warnState = warning('off','all');

lb = -10*ones(1, na+nc);
ub =  10*ones(1, na+nc);

if isempty(gcp('nocreate'))
    % start pool without printing by capturing output
    evalc('parpool'); 
end

obj = @(p) arma_loss(p, z, na, nc);
[popt, Jmin] = fmincon(obj, p0, [], [], [], [], lb, ub, [], opts);

% Reconstruct final ARMA model
Aopt = [1, popt(1:na)];
Copt = [1, popt(na+1:na+nc)];
model1 = idpoly(Aopt, [], Copt);

% fprintf('Optimized loss = %.6g\n', Jmin);

% Restore warning state
warning(warnState);

end