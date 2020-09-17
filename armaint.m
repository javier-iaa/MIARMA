function [interp, go] = armaint(seg1, seg2, ord, N2)
% function [interp,go] = armaint(seg1, seg2, ord, N2) interpolates N2
% data points between the segments seg1 and seg2 using ARMA models.
% To generate the output segment interp a triangular weight is used for
% both segments.
% Inputs:       seg1 - left data segment
%               seg2 - right data segment
%               ord - ARMA (p,q) orders
%               N2 - length of the gap
% Outputs:      interp - interpolated segment
%               go - true when the interpolation works and false otherwise
% Version: 1.3.9
% Changes from the last version:
% - Removed version control
% - Minor improvements.
%
%  Calls: sigma_clip.m
%  Author(s): Javier Pascual-Granado
%  Date: 17/09/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

go = true;
interp = NaN;
msg = [];

% Coefficient used to detect if the extrapolation explodes
fac_sig = 5;

%% Algorithm properties
% Zstability: Specifies the maximum distance of all poles 
% from the origin to test stability of discrete-time models. 
% A model is considered stable if all poles are within the 
% distance Zstability from the origin. Default is 1+sqrt(eps)
stab = 1+sqrt(eps);
myalg.Focus = 'Prediction'; % Prediction, Simulation or Stability
myalg.MaxIter = 20;
% myalg.Tolerance = 0.0100; % default value
myalg.Tolerance = 0.05;
myalg.LimitError = 0; % Specifies when to adjust the weight of large errors 
% from quadratic to linear. Default value is 0. Errors larger than LimitError 
% times the estimated standard deviation have a linear weight in the criteria.
myalg.MaxSize = 'Auto'; % data is split into segments where each contains fewer than MaxSize elements.
% myalg.SearchMethod = 'Auto';
myalg.SearchMethod = 'lsqnonlin'; % options are gn, gna, lm, Auto, lsqnonlin
myalg.Criterion = 'Det';    % Det or Trace
myalg.Weighting = 1;
myalg.FixedParameter = [];
myalg.Display = 'Off';
myalg.N4Weight = 'Auto';
myalg.N4Horizon = 'Auto';
% myalg.InitialState = 'Backcast';
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
    ts = 1;
    data2 = iddata(flipud(seg2n(2:end)),[],ts);
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
        ts = 1;
        data2 = iddata(flipud(seg2n(2:end)),[],ts);
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
    model1 = armax(seg1n,ord,'alg',myalg);
catch E
    msg = getReport(E);
    go = false;
    return
end
ts = 1;
data1 = iddata(seg1n(1:end-1),[],ts);
yfor = pred(model1,data1,N2+1,'e');
%     yfor = pred(model1,data1,N2,'e');
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
ts = 1;
data2 = iddata(flipud(seg2n(2:end)),[],ts);
yback = pred(model2,data2,N2+1,'e');
%     yback = pred(model2,data2,N2,'e');
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

if sigd_yf > fac_sig*sigd_s1 && sigd_yf > fac_sig*sigd_s2
    
    seg1 = sigma_clip(seg1, 2);
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
    ts = 1;
    data1 = iddata(seg1n(1:end-1),[],ts);
    yfor = pred(model1,data1,N2+1,'e');
%         yfor = pred(model1,data1,N2,'e');
    yfor = yfor.y;
          
    % Estimation of residuals
    e1 = resid(model1, data1);
    rstd1 = std(e1.OutputData);
    r1 = rstd1*randn(size(yfor));
    yfor = ( (yfor + r1).*sig_s1 + mean(seg1) );
    
    sigd_yf = std(diff(yfor));
    sigd_s1 = std(diff(seg1));
    
    if sigd_yf > fac_sig*sigd_s1 && sigd_yf > fac_sig*sigd_s2
%         yfor = yfor.*(sig_s1/sig_yf);
       interp = NaN(1,N2);
       go = false;
       return
    end
end

% Validation with sigma clipping
% This code can produce bugs and should be used with care and only when
% the extrapolations are unstable.

sigd_yb = std(diff(yback));

if sigd_yb > fac_sig*sigd_s2 && sigd_yb > fac_sig*sigd_s1
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
    ts = 1;
    data2 = iddata(flipud(seg2n(2:end)),[],ts);
    yback = pred(model2,data2,N2+1,'e');
%         yback = pred(model2,data2,N2,'e');
    yback = yback.y;
    yback = flipud(yback);
    
    % Estimation of residuals
    e2 = resid(model2, data2);
    rstd2 = std(e2.OutputData);
    r2 = rstd2*randn(size(yback));
    yback = ( (yback + r2).*sig_s2 + mean(seg2) );
    
    sigd_yb = std(diff(yback));
    sigd_s2 = std(diff(seg2));
    
    if sigd_yb > fac_sig*sigd_s2 && sigd_yb > fac_sig*sigd_s1
%         yback = yback.*(sig_s2/sig_yb);
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
