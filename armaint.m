function [interp, go] = armaint(seg1, seg2, ord, N2)
% function [interp,go] = armaint(seg1, seg2, ord, N2) interpolates N2
% data points between the segments seg1 and seg2 using ARMA models.
% To generate the output segment interp a triangular weight is used for
% both segments.
% Inputs:   seg1 - left data segment
%               seg2 - right data segment
%               ord - ARMA (p,q) orders
%               N2 - length of the gap
% Outputs:      interp - interpolated segment
%               go - true when the interpolation works and false otherwise
% Version: 1.4.1
% Changes from the last version:
% - Implemented a nan check for yfor and yback
% - Crash test avoid 'compare' failing for long predictions.
% - Conditions for sigma clipping are explicit.
% - Added alternative options for armax in validation
%
%  Calls: sigma_clip.m
%  Author(s): Javier Pascual-Granado
%  Date: 12/01/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

go = true;
interp = NaN;
msg = [];

% Coefficient used to detect if the extrapolation explodes
fac_sig = 5;

% Limit in goodness of fitting of models to continue with extrap
lim_gf = 80;

%% Algorithm properties
% Zstability: Specifies the maximum distance of all poles 
% from the origin to test stability of discrete-time models. 
% A model is considered stable if all poles are within the 
% distance Zstability from the origin. Default is 1+sqrt(eps)
stab = 1+sqrt(eps);
myalg.Focus = 'Stability'; % Prediction, Simulation or Stability
myalg.MaxIter = 20;
% myalg.Tolerance = 0.0100; % default value
myalg.Tolerance = 0.05;
myalg.LimitError = 0; % Specifies when to adjust the weight of large errors 
% from quadratic to linear. Default value is 0. Errors larger than LimitError 
% times the estimated standard deviation have a linear weight in the criteria.
myalg.MaxSize = 'Auto'; % data is split into segments where each contains fewer than MaxSize elements.
% myalg.SearchMethod = 'Auto';
myalg.SearchMethod = 'lm'; % options are gn, gna, lm, Auto, lsqnonlin
myalg.Criterion = 'Det';    % Det or Trace
myalg.Weighting = 1;
myalg.FixedParameter = [];
myalg.Display = 'Off';
myalg.N4Weight = 'Auto';
myalg.N4Horizon = 'Auto';
myalg.InitialState = 'Auto';
myalg.Advanced.Search.GnPinvConst = 10000;
myalg.Advanced.Search.InitGnaTol = 1.0e-04;
myalg.Advanced.Search.LmStep = 2;
myalg.Advanced.Search.StepReduction = 2;
myalg.Advanced.Search.MaxBisections = 25;
myalg.Advanced.Search.LmStartValue= 1.0e-03;
myalg.Advanced.Search.RelImprovement = 0;
myalg.Advanced.Threshold.Zstability = stab; % Specifies the maximum distance 
% of all poles from the origin to test stability of discrete-time models.
% default is 1+sqrt(eps)
myalg.Advanced.Threshold.Sstability = 0; % Specifies the location of the 
% rightmost pole to test the stability of continuous-time models. Default
% is zero.
myalg.Advanced.Threshold.AutoInitialState = 1.05; %  Specifies when to automatically
% estimate the initial state. When InitialState = 'Auto', the initial state is 
% estimated when the ratio of the prediction-error norm with a zero initial state 
% to the norm with an estimated initial state exceeds AutoInitialState. 
% Default is 1.05.

%% Preparing data
% Stationarity is assumed
sigma = std([seg1; seg2]);
if isinf(sigma)==1
    fprintf('Infinite deviation error\n');
    interp = NaN(1,N2);
    go = false;
    return;
end

% Change row vectors into column vectors
if ~isempty(find(isnan(seg2),1))
    [fil,~] = size(seg1);
    if fil==1
        seg1 = seg1';
    end
elseif ~isempty(find(isnan(seg1),1))
    [fil,~] = size(seg2);
    if fil==1
        seg2 = seg2';
    end
end

% Weights
wp = 1/(N2+1);
w = (wp:wp:(1-wp));

%% Forward predictor: ARMA approach using an iterative algorithm
if ~isempty(find(isnan(seg2),1))
    sig_s1 = std(seg1);
    
    % Normalization
    seg1n = (seg1-mean(seg1))./sig_s1;

    % Calculate ARMA model and obtain the coeff. for the left segment
    try
        model1 = armax(seg1n,ord,'alg',myalg);
    catch E
        go = false;
        msg = getReport(E);
        return
    end
    ts = 1;
    data1 = iddata(seg1n(1:end-1),[],ts);
    yfor = pred(model1,data1,N2+1,'e');
    yfor = yfor.y;
        
    % Estimation of residuals
    e1 = resid(model1, data1);
    rstd1 = std(e1.OutputData);
    r1 = rstd1*randn(size(yfor));
    yfor = ( (yfor + r1).*sig_s1 + mean(seg1) );

% Validation with sigma clipping
% This code can produce bugs and should be used with care and only when
% the extrapolations are unstable.
    
    sig_yf = std(yfor);

    if sig_yf > fac_sig*sig_s1
        seg1 = sigma_clip(seg1,2);
        seg1n = (seg1-mean(seg1))./std(seg1);
        sig_s1 = std(seg1);

        % Calculate ARMA model and obtain the coeff. for the left segment
        try
            model1 = armax(seg1n,ord,'alg',myalg);
        catch E
            msg = getReport(E);
            go = false;
            return
        end
        data1 = iddata(seg1n(1:end-1),[],ts);
        yfor = pred(model1,data1,N2+1,'e');
        yfor = yfor.y;
            
        % Estimation of residuals
        e1 = resid(model1, data1);
        rstd1 = std(e1.OutputData);
        r1 = rstd1*randn(size(yfor));
        yfor = ( (yfor + r1).*sig_s1 + mean(seg1) );

        sig_yf = std(yfor);

        if sig_yf > fac_sig*sig_s1
            interp = NaN(1,N2);
            go = false;
            return
        end
    end

    % Interpolated data
    interp = yfor(2:end);
    return
end

%% Backward predictor: ARMA approach using an iterative algorithm
if ~isempty(find(isnan(seg1),1))
    sig_s2 = std(seg2);
    
    % Normalization
    seg2n = (seg2-mean(seg2))./sig_s2;

    % Calculate ARMA model and obtain the coeff. for the right segment
    try
        model2 = armax(flipud(seg2n), ord,'alg',myalg);
    catch E
        msg = getReport(E);
        go = false;
        return
    end
    data2 = iddata( flipud(seg2n(2:end)), [] );
    yback = pred(model2,data2,N2+1,'e');
    yback = yback.y;
    yback = flipud(yback);

    % Estimation of residuals
    e2 = resid(model2, data2);
    rstd2 = std(e2.OutputData);
    r2 = rstd2*randn(size(yback));
    yback = ( (yback + r2).*sig_s2 + mean(seg2) );

% Validation with sigma clipping
% This code can produce bugs and should be used with care and only when
% the extrapolations are unstable.
    
    sig_yb = std(yback);

    if sig_yb > fac_sig*sig_s2
        seg2 = sigma_clip(seg2,2);
        seg2n = (seg2-mean(seg2))./std(seg2);
        sig_s2 = std(seg2);

        % Calculate ARMA model and obtain the coeff. for the right segment
        try
            model2 = armax(flipud(seg2n), ord,'alg',myalg);
        catch E
            msg = getReport(E);
            go = false;
            return
        end
        data2 = iddata( flipud(seg2n(2:end)) ,[] );
        yback = pred(model2,data2,N2+1,'e');
        yback = yback.y;
        yback = flipud(yback);
        
        % Estimation of residuals
        e2 = resid(model2, data2);
        rstd2 = std(e2.OutputData);
        r2 = rstd2*randn(size(yback));
        yback = ( (yback + r2).*sig_s2 + mean(seg2) );
        
        sig_yb = std(yback);

        if sig_yb > fac_sig*sig_s2
            interp = NaN(1,N2);
            go = false;
            return
        end
    end

    % Interpolated data
    interp = yback(1:end-1);
    return
end

%% Forward-Backward predictor: ARMA approach using an iterative algorithm

% Alternative options
myalg_alt = myalg;
myalg_alt.Focus = 'Prediction';

% Normalization
sig_s1 = std(seg1);
sig_s2 = std(seg2);
seg1n = (seg1-mean(seg1))./sig_s1;
seg2n = (seg2-mean(seg2))./sig_s2;

% Calculate ARMA model, obtain the coefficients and validate 
% forecasted data through a 2-sigma criterion

% Forward extrapolation
% Calculate ARMA model and obtain the coeff. for the left segment
try
    model1 = armax(seg1n, ord, 'alg', myalg);
catch E
    msg = getReport(E);
    go = false;
    return
end

data1 = iddata( seg1n(1:end-1), []);
yfor = pred( model1, data1, N2+1, 'e' );
yfor = yfor.y;

% Backward extrapolation

% Calculate ARMA model and obtain the coeff. for the right segment
try
    model2 = armax(flipud(seg2n), ord,'alg',myalg);
catch E
    msg = getReport(E);
    go = false;
    return
end

data2 = iddata( flipud(seg2n(2:end)), [] );
yback = pred(model2,data2,N2+1,'e');
yback = yback.y;
yback = flipud(yback);

% Estimation of residuals
e1 = resid(model1, data1);
rstd1 = std(e1.OutputData);
r1 = rstd1*randn(size(yfor));
yfor = ( (yfor + r1).*sig_s1 + mean(seg1) );
e2 = resid(model2, data2);
rstd2 = std(e2.OutputData);
r2 = rstd2*randn(size(yback));
yback = ( (yback + r2).*sig_s2 + mean(seg2) );

% Validation with sigma clipping
% This code should be used with care when the extrapolations are unstable.

% The sigma of the differences shows better when the extrapolation explodes
sigd_yf = std(diff(yfor));
sigd_s1 = std(diff(seg1));
sigd_s2 = std(diff(seg2));

% It seems utcompare function implemented in Matlab breaks for longer predictions
lseg = length(seg1);
qL = 0.8; % limit for prediction horizon
predh = floor( qL*lseg );

% Goodness of fitting for left model
if N2>lseg,
    [~,g1,~] = compare(data1, model1, predh);
else
    try
        [~,g1,~] = compare(data1, model1, N2+1);
    catch E
        % This is for debugging purposes and will be eliminated in later versions
        fprintf('\nN2 = %d, d = %d, L = %d ... \n', N2, sum(ord), lseg );
        go = false;
        return
    end
%     [~,g1,~] = compare(data1, model1, N2+1);
end    


% Conditions to enable sigma clipping
cfcom = g1<lim_gf;
cf1sig = sigd_yf > fac_sig*sigd_s1;
cf2sig = sigd_yf > fac_sig*sigd_s2;
cfsig = cf1sig && cf2sig && cfcom;
connanf = isnan(sigd_yf);

if  cfsig || connanf
    seg1 = sigma_clip( seg1 );
    sig_s1 = std( seg1 );
    seg1n = (seg1-mean(seg1))./sig_s1;
    
    % Calculate ARMA model and obtain the coeff. for the left segment
    try
        model1 = armax(seg1n, ord, 'alg', myalg_alt);
    catch E
        msg = getReport(E);
        go = false;
        return
    end
    data1 = iddata( seg1n(1:end-1), [] );
    yfor = pred( model1, data1, N2+1, 'e' );
    yfor = yfor.y;
          
    % Estimation of residuals
    e1 = resid(model1, data1);
    rstd1 = std(e1.OutputData);
    r1 = rstd1*randn(size(yfor));
    yfor = ( (yfor + r1).*sig_s1 + mean(seg1) );
    
    sigd_yf = std(diff(yfor));
    sigd_s1 = std(diff(seg1));
    if N2>lseg,
        [~,g1,~] = compare(data1, model1, predh);
    else
        [~,g1,~] = compare(data1, model1, N2+1);
    end 
    
    % Conditions for sigma clipping
    cfcom = g1<lim_gf;
    cf1sig = sigd_yf > fac_sig*sigd_s1;
    cf2sig = sigd_yf > fac_sig*sigd_s2;
    cfsig = cf1sig && cf2sig && cfcom;
    connanf = isnan(sigd_yf);

    if  cfsig || connanf
       interp = NaN(1,N2);
       go = false;
       return
    end
end

% Validation with sigma clipping
% This code can produce bugs and should be used with care and only when
% the extrapolations are unstable.

sigd_yb = std(diff(yback));
if N2>lseg,
    [~,g2,~] = compare(data2, model2, predh);
else
    [~,g2,~] = compare(data2, model2, N2+1);
end

% Conditions to enable sigma clipping
cbcom = g2<lim_gf;
cb1sig = sigd_yb > fac_sig*sigd_s1;
cb2sig = sigd_yb > fac_sig*sigd_s2;
cbsig = cb1sig && cb2sig && cbcom;
cbnan = isnan(sigd_yb);

if  cbsig || cbnan
    seg2 = sigma_clip( seg2 );
    sig_s2 = std(seg2);
    seg2n = (seg2-mean(seg2))./sig_s2;
    
    % Calculate ARMA model and obtain the coeff. for the right segment
    try
        model2 = armax( flipud(seg2n), ord, 'alg', myalg_alt );
    catch E
        msg = getReport(E);
        go = false;
        return
    end
    data2 = iddata( flipud(seg2n(2:end)), [] );
    yback = pred( model2, data2, N2+1, 'e');
    yback = yback.y;
    yback = flipud(yback);
    
    % Estimation of residuals
    e2 = resid(model2, data2);
    rstd2 = std(e2.OutputData);
    r2 = rstd2*randn(size(yback));
    yback = ( (yback + r2).*sig_s2 + mean(seg2) );
    
    sigd_yb = std(diff(yback));
    sigd_s2 = std(diff(seg2));
    if N2>lseg,
        [~,g2,~] = compare(data2, model2, predh);
    else
        [~,g2,~] = compare(data2, model2, N2+1);
    end 

    % Conditions to enable sigma clipping
    cbcom = g2<lim_gf;
    cb1sig = sigd_yb > fac_sig*sigd_s1;
    cb2sig = sigd_yb > fac_sig*sigd_s2;
    cbsig = cb1sig && cb2sig && cbcom;
    cbnan = isnan(sigd_yb);

    if  cbsig || cbnan
        interp = NaN(1,N2);
        go = false;
        return
    end
end

% Interpolated data
if isempty(msg)
    interp = (1-w').*yfor(2:end) + w'.*yback(1:end-1);
%     interp = (1-w').*yfor + w'.*yback;
  % interp = interp + rstd*randn(size(interp));
end

end
