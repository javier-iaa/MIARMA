function [interp,go] = armaint(seg1,seg3,ord,N2)
% function [interp,go] = armaint(seg1,seg3,ord,N2) interpolates N2
% data points between the segments seg1 and seg3 using ARMA models.
% To generate the output segment interp a triangular weight is used for
% both segments.
% Inputs:       seg1 - left data segment
%               seg3 - right data segment
%               ord - ARMA (p,q) orders
%               N2 - length of the gap
% Outputs:      interp - interpolated segment
%               go - true when the interpolation works and false otherwise

% Version: 1.3.4
% Changes from the last version:
% Estimation of rstd through fitting residuals.
%
%  Calls: sigma_clip.m
%  Author(s): Javier Pascual-Granado
%  $Date: 16/10/2019$
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

go = true;

% Std factor for sigma clipping
fac_sig = 3;

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

%% Version control
% Three different versions are available: old, new and default (nothing is
% done).
rel = version;
rel = str2double(rel(end-5:end-2));
modo = 'default';

if rel < 2012
    modo = 'old';
elseif rel >= 2015
    modo = 'new';
end

%% Preparing data
% Stationarity is assumed
sigma = std([seg1; seg3]);
if isinf(sigma)==1
    fprintf('Infinite deviation error\n');
    interp = NaN(1,N2);
    go = false;
    return;
end

% Change row vectors into column vectors
if ~isempty(find(isnan(seg3),1))
    [fil,~] = size(seg1);
    if fil==1
        seg1 = seg1';
    end
elseif ~isempty(find(isnan(seg1),1))
    [fil,~] = size(seg3);
    if fil==1
        seg3 = seg3';
    end
end

% Weights
wp = 1/(N2+1);
w = (wp:wp:(1-wp));

%% Forward predictor: ARMA approach using an iterative algorithm
if ~isempty(find(isnan(seg3),1))
    sig_s1 = std(seg1);
    
    % Normalization
    seg1n = (seg1-mean(seg1))./sig_s1;

    if strcmp(modo, 'old')
        % Calculate ARMA model and obtain the coeff. for the left segment
        model1 = armax(seg1n,ord,'alg',myalg);
        ts = 1;
        data1 = iddata(seg1n(1:end-1),[],ts);
        yfor = pred(model1,data1,N2+1,'e');
        yfor = yfor.y;
        
    elseif strcmp(modo, 'new')
        % Set the options
        opt = armaxOptions;
        opt.Advanced.ErrorThreshold = 1.6; % put zero here if no outliers
        opt.Display = 'off';
        % Then fit the model
        model1 = armax(seg1n, ord, opt);
        yfor = forecast(model1,seg1n(1:end-1),N2+1);
        
    else
       % Calculate ARMA model and obtain the coeff. for the left segment
        model1 = armax(seg1n, ord,'alg',myalg);
        yfor = forecast(model1, seg1n(1:end-1), N2+1);
    end

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

        if strcmp(modo,'old')
            % Calculate ARMA model and obtain the coeff. for the left segment
            model1 = armax(seg1n,ord,'alg',myalg);
            data1 = iddata(seg1n(1:end-1),[],ts);
            yfor = pred(model1,data1,N2+1,'e');
            yfor = yfor.y;
            
        elseif strcmp(modo,'new')
            % Set the options
            opt = armaxOptions;
            opt.Advanced.ErrorThreshold = 1.6; % put zero here if no outliers
            opt.Display = 'off';
            % Then fit the model
            model1 = armax(seg1n, ord, opt);
            yfor = forecast(model1,seg1n(1:end-1),N2+1);
            
        else
            % Calculate ARMA model and obtain the coeff. for the left segment
            model1 = armax(seg1n,ord,'alg',myalg);
            yfor = forecast(model1,seg1n(1:end-1),N2+1); 
        end

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
    interp = interp + rstd*randn(size(interp)); 
    return
end
%% Backward predictor: ARMA approach using an iterative algorithm
if ~isempty(find(isnan(seg1),1))
    sig_s3 = std(seg3);
    
    % Normalization
    seg3n = (seg3-mean(seg3))./sig_s3;

    if strcmp(modo, 'old')
        % Calculate ARMA model and obtain the coeff. for the right segment
        model2 = armax(flipud(seg3n), ord,'alg',myalg);
        ts = 1;
        data2 = iddata(flipud(seg3n(2:end)),[],ts);
        yback = pred(model2,data2,N2+1,'e');
        yback = yback.y;
        
    elseif strcmp(modo, 'new')
        % Set the options
        opt = armaxOptions;
        opt.Advanced.ErrorThreshold = 1.6; % put zero here if no outliers
        opt.Display = 'off';
        % Then fit the model
        model2 = armax(flipud(seg3n), ord, opt);
        yback = forecast(model2, flipud(seg3n(2:end)), N2+1);

    else
        
        % Calculate ARMA model and obtain the coeffs for the right segment
        model2 = armax(flipud(seg3n),ord,'alg',myalg);
        yback = forecast(model2,flipud(seg3n(2:end)),N2+1);
    end

    yback = flipud(yback);

    % Estimation of residuals
    e2 = resid(model2, data2);
    rstd2 = std(e2.OutputData);
    r2 = rstd2*randn(size(yback));
    yback = ( (yback + r2).*sig_s3 + mean(seg3) );

% Validation with sigma clipping
% This code can produce bugs and should be used with care and only when
% the extrapolations are unstable.
    
    sig_yb = std(yback);

    if sig_yb > fac_sig*sig_s3
        seg3 = sigma_clip(seg3,2);
        seg3n = (seg3-mean(seg3))./std(seg3);
        sig_s3 = std(seg3);

        if strcmp(modo, 'old')
            % Calculate ARMA model and obtain the coeff. for the right segment
            model2 = armax(flipud(seg3n), ord,'alg',myalg);
            ts = 1;
            data2 = iddata(flipud(seg3n(2:end)),[],ts);
            yback = pred(model2,data2,N2+1,'e');
            yback = yback.y;

        elseif strcmp(modo, 'new')
            % Set the options
            opt = armaxOptions;
            opt.Advanced.ErrorThreshold = 1.6; % put zero here if no outliers
            opt.Display = 'off';
            % Then fit the model
            model2 = armax(flipud(seg3n), ord, opt);
            yback = forecast(model2, flipud(seg3n(2:end)), N2+1);

        else
            % Calculate ARMA model and obtain the coeffs for the right segment
            model2 = armax(flipud(seg3n),ord,'alg',myalg);
            yback = forecast(model2,flipud(seg3n(2:end)),N2+1);
        end

        yback = flipud(yback);
        
        % Estimation of residuals
        e2 = resid(model2, data2);
        rstd2 = std(e2.OutputData);
        r2 = rstd2*randn(size(yback));
        yback = ( (yback + r2).*sig_s3 + mean(seg3) );
        
        sig_yb = std(yback);

        if sig_yb > fac_sig*sig_s3
            interp = NaN(1,N2);
            go = false;
            return
        end
    end

    % Interpolated data
    interp = yback(1:end-1);
    interp = interp + rstd*randn(size(interp));
    return
end

%% Forward-Backward predictor: ARMA approach using an iterative algorithm
% Normalization
sig_s1 = std(seg1);
sig_s3 = std(seg3);
seg1n = (seg1-mean(seg1))./sig_s1;
seg3n = (seg3-mean(seg3))./sig_s3;

% Calculate ARMA model, obtain the coefficients and validate 
% forecasted data through a 2-sigma criterion

% Forward extrapolation
if strcmp(modo, 'old')
    % Calculate ARMA model and obtain the coeff. for the left segment
    model1 = armax(seg1n,ord,'alg',myalg);
    ts = 1;
    data1 = iddata(seg1n(1:end-1),[],ts);
    yfor = pred(model1,data1,N2+1,'e');
    yfor = yfor.y;

elseif strcmp(modo, 'new')
    % Set the options
    opt = armaxOptions;
    opt.Advanced.ErrorThreshold = 1.6; % put zero here if no outliers
    opt.Display = 'off';
    % Then fit the model
    model1 = armax(seg1n, ord, opt);
    yfor = forecast(model1,seg1n(1:end-1),N2+1);

else
   % Calculate ARMA model and obtain the coeff. for the left segment
    model1 = armax(seg1n, ord,'alg',myalg);
    yfor = forecast(model1, seg1n(1:end-1), N2+1);
end

% Backward extrapolation
if strcmp(modo, 'old')
    % Calculate ARMA model and obtain the coeff. for the right segment
    model2 = armax(flipud(seg3n), ord,'alg',myalg);
    ts = 1;
    data2 = iddata(flipud(seg3n(2:end)),[],ts);
    yback = pred(model2,data2,N2+1,'e');
    yback = yback.y;

elseif strcmp(modo, 'new')
    % Set the options
    opt = armaxOptions;
    opt.Advanced.ErrorThreshold = 1.6; % put zero here if no outliers
    opt.Display = 'off';
    % Then fit the model
    model2 = armax(flipud(seg3n), ord, opt);
    yback = forecast(model2, flipud(seg3n(2:end)), N2+1);

else
    % Calculate ARMA model and obtain the coeffs for the right segment
    model2 = armax(flipud(seg3n),ord,'alg',myalg);
    yback = forecast(model2,flipud(seg3n(2:end)),N2+1);
end

yback = flipud(yback);

% Estimation of residuals
e1 = resid(model1, data1);
rstd1 = std(e1.OutputData);
r1 = rstd1*randn(size(yfor));
yfor = ( (yfor + r1).*sig_s1 + mean(seg1) );
e2 = resid(model2, data2);
rstd2 = std(e2.OutputData);
r2 = rstd2*randn(size(yback));
yback = ( (yback + r2).*sig_s3 + mean(seg3) );

% Validation with sigma clipping
% This code should be used with care when the extrapolations are unstable.

sig_yf = std(yfor);

if sig_yf > fac_sig*sig_s1 && sig_yf > fac_sig*sig_s3
    seg1 = sigma_clip(seg1, 2);
    seg1n = (seg1-mean(seg1))./std(seg1);
    sig_s1 = std(seg1);
    
    if strcmp(modo, 'old')
        % Calculate ARMA model and obtain the coeff. for the left segment
        model1 = armax(seg1n,ord,'alg',myalg);
        ts = 1;
        data1 = iddata(seg1n(1:end-1),[],ts);
        yfor = pred(model1,data1,N2+1,'e');
        yfor = yfor.y;
        
    elseif strcmp(modo, 'new')
        % Set the options
        opt = armaxOptions;
        opt.Advanced.ErrorThreshold = 1.6; % put zero here if no outliers
        opt.Display = 'off';
        % Then fit the model
        model1 = armax(seg1n, ord, opt);
        yfor = forecast(model1,seg1n(1:end-1),N2+1);
        
    else
       % Calculate ARMA model and obtain the coeff. for the left segment
        model1 = armax(seg1n, ord,'alg',myalg);
        yfor = forecast(model1, seg1n(1:end-1), N2+1);
    end
    
    % Estimation of residuals
    e1 = resid(model1, data1);
    rstd1 = std(e1.OutputData);
    r1 = rstd1*randn(size(yfor));
    yfor = ( (yfor + r1).*sig_s1 + mean(seg1) );
    
    sig_yf = std(yfor);
    
    if sig_yf > fac_sig*sig_s1 && sig_yf > fac_sig*sig_s3
        yfor = yfor.*(sig_s1/sig_yf);
%        interp = NaN(1,N2);
%        go = false;
%        return
    end
end

% Validation with sigma clipping
% This code can produce bugs and should be used with care and only when
% the extrapolations are unstable.

sig_yb = std(yback);

if sig_yb > fac_sig*sig_s3 && sig_yb > fac_sig*sig_s1
    seg3 = sigma_clip(seg3,2);
    seg3n = (seg3-mean(seg3))./std(seg3);
    sig_s3 = std(seg3);
    
    if strcmp(modo, 'old')
        % Calculate ARMA model and obtain the coeff. for the right segment
        model2 = armax(flipud(seg3n), ord,'alg',myalg);
        ts = 1;
        data2 = iddata(flipud(seg3n(2:end)),[],ts);
        yback = pred(model2,data2,N2+1,'e');
        yback = yback.y;
        
    elseif strcmp(modo, 'new')
        % Set the options
        opt = armaxOptions;
        opt.Advanced.ErrorThreshold = 1.6; % put zero here if no outliers
        opt.Display = 'off';
        % Then fit the model
        model2 = armax(flipud(seg3n), ord, opt);
        yback = forecast(model2, flipud(seg3n(2:end)), N2+1);

    else   
        % Calculate ARMA model and obtain the coeffs for the right segment
        model2 = armax(flipud(seg3n),ord,'alg',myalg);
        yback = forecast(model2,flipud(seg3n(2:end)),N2+1);
    end
    yback = flipud(yback);
    
    % Estimation of residuals
    e2 = resid(model2, data2);
    rstd2 = std(e2.OutputData);
    r2 = rstd2*randn(size(yback));
    yback = ( (yback + r2).*sig_s3 + mean(seg3) );
    
    sig_yb = std(yback);
    
    if sig_yb > fac_sig*sig_s3 && sig_yb > fac_sig*sig_s1
        yback = yback.*(sig_s3/sig_yb);
%         interp = NaN(1,N2);
%         go = false;
%         return
    end
end

% Interpolated data
interp = (1-w').*yfor(2:end) + w'.*yback(1:end-1);
% interp = interp + rstd*randn(size(interp));

end