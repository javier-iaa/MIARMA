function [interp, go, info] = armaint(seg1, seg2, ord, N2, varargin)
% function [interp,go, info] = armaint(seg1, seg2, ord, N2, varargin) 
% interpolates N2 data points between the segments seg1 and seg2 using ARMA
%  models.
% To generate the output segment interp a triangular weight is used for
% both segments.
% Inputs:       seg1 - left data segment
%               seg2 - right data segment
%               ord - ARMA (p,q) orders
%               N2 - length of the gap
%               
% Optional inputs:
%               mem - can be used to limit memory use. Pass the flag 'mem'
%                followed by the number in Gb.
%               debug - is used for debug/test purposes. Pass the flag
%                'debug' to activate it.
%
% Outputs:      interp - interpolated segment
%               go - true when the interpolation works and false otherwise
%               info - optional output structure containing information of
%               the performance of armax algorithms.
%
% Version: 1.4.7 - R2024
%
% Changes from the last version:
% - Reverted armax_par to armax.
%
%  Calls: sigma_clip.m, algoprop.m
%  Author(s): Javier Pascual-Granado
%  Date: 30/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This flag load the customised algorithm options included in algoprop,
% otherwise, default options are used. I leave this activated by
% construction but this should change at some point.
myalg_flag = true;

if myalg_flag
    myalg = algoprop();
end

% Set the debug flag for testing purposes
debug_flag = find(strcmp(varargin,'debug'), 1 );
if ~isempty(debug_flag)
    debug = true;
else
    debug = false;
end

% Default value for mem in Gb
mem_flag = find(strcmp(varargin,'mem'), 1 );
if ~isempty(mem_flag)
    mem = varargin{mem_flag+1};
else
    mem = 16;
end

go = true;
interp = NaN;
msg = [];

% Coefficient used to detect if the extrapolation explodes
fac_sig = 5;

% Limit in goodness of fitting of models to continue with extrap
lim_gf = 80;

%% Forecast options
opt = forecastOptions('InitialCondition', 'e');

%% Preparing data
% Stationarity is assumed
sigma = std([seg1; seg2]);
if isinf(sigma)
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

po = ord(1);

% This is approx. the segment limit to avoid a memory overflow
lim_segsize = floor( mem/(8*(po^3 + 4*po^2 + po)/1024^3) - 10 );
lim_segsize = lim_segsize - 10; % small adjust to avoid overflow

if length(seg1) > lim_segsize
    seg1 = tail(seg1, lim_segsize);
end

if length(seg2) > lim_segsize
    seg2 = head(seg2, lim_segsize);
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
        % model1 = armax_par(seg1n,ord);
        model1 = armax(seg1n,ord);
    catch E
        go = false;
        msg = getReport(E);
        return
    end
    ts = 1;
    data1 = iddata(seg1n(1:end-1),[],ts);
    yfor = forecast(model1,data1,N2+1, opt);
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
            % model1 = armax_par(seg1n,ord);
            model1 = armax(seg1n,ord);
        catch E
            msg = getReport(E);
            go = false;
            return
        end
        data1 = iddata(seg1n(1:end-1),[],ts);
        yfor = forecast(model1,data1,N2+1,opt);
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
        % model2 = armax_par(seg2n,ord);
        model2 = armax(seg2n,ord);
    catch E
        msg = getReport(E);
        go = false;
        return
    end
    data2 = iddata( flipud(seg2n(2:end)), [] );
    yback = forecast(model2,data2,N2+1,opt);
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
            % model2 = armax_par(seg2n,ord);
            model2 = armax(seg2n,ord);
        catch E
            msg = getReport(E);
            go = false;
            return
        end
        data2 = iddata( flipud(seg2n(2:end)) ,[] );
        yback = forecast(model2,data2,N2+1,opt);
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

% Alternative options (with focus on prediction)
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

% This is used for debugging/test purposes. In the future it will be
% removed and only one of the armax algorithms will be preserved. Then, it
% should be added to the forward and backward extrapolation sections too.
if debug
    % Input info
    info = struct();
    info.numel = numel(seg1n);
    info.std = sig_s1;
    info.ordp = ord(1);
    info.ordq = ord(2);
    
    % =========================================================
    % 1. armax_par(seg1n,ord)
    % =========================================================
    tic;
    try
        [model1, info.armax_par_loss] = armax_par(seg1n,ord);
        info.armax_par_time = toc;
        info.armax_par_ok   = true;
    catch E
        info.armax_par_loss = nan;
        info.armax_par_time = toc;
        info.armax_par_ok   = false;
        info.armax_par_error = getReport(E,'basic');
        model1 = [];
    end

    % =========================================================
    % 2. armax(seg1n,ord)
    % =========================================================
    tic;
    try
        model12 = armax(seg1n,ord);
        info.armax_loss = model12.Report.Fit.LossFcn;
        info.armax_time = toc;
        info.armax_ok   = true;
    catch E
        info.armax_loss = nan;
        info.armax_time = toc;
        info.armax_ok   = false;
        info.armax_error = getReport(E,'basic');
    end

    % =========================================================
    % 3. armax(seg1n,ord,'alg',myalg)
    % =========================================================
    tic;
    try
        model13 = armax(seg1n,ord,'alg',myalg);
        info.armax_alg_loss = model13.Report.Fit.LossFcn;
        info.armax_alg_time = toc;
        info.armax_alg_ok   = true;
    catch E
        info.armax_alg_loss = nan;
        info.armax_alg_time = toc;
        info.armax_alg_ok   = false;
        info.armax_alg_error = getReport(E,'basic');
    end

else
    % Normal execution
    try
        % model1 = armax_par(seg1n,ord);
        model1 = armax(seg1n,ord);
    catch E
        msg = getReport(E);
        return
    end
end

data1 = iddata( seg1n(1:end-1), []);
yfor = forecast( model1, data1, N2+1, opt );
yfor = yfor.y;

% Backward extrapolation

% Calculate ARMA model and obtain the coeff. for the right segment
try
    % model2 = armax_par(seg2n,ord);
    model2 = armax(seg2n,ord);
catch E
    msg = getReport(E);
    go = false;
    return
end

data2 = iddata( flipud(seg2n(2:end)), [] );
yback = forecast(model2,data2,N2+1,opt);
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
if N2>lseg
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
        % model1 = armax_par(seg1n,ord);
        model1 = armax(seg1n,ord);
    catch E
        msg = getReport(E);
        go = false;
        return
    end
    data1 = iddata( seg1n(1:end-1), [] );
    yfor = forecast( model1, data1, N2+1, opt );
    yfor = yfor.y;
          
    % Estimation of residuals
    e1 = resid(model1, data1);
    rstd1 = std(e1.OutputData);
    r1 = rstd1*randn(size(yfor));
    yfor = ( (yfor + r1).*sig_s1 + mean(seg1) );
    
    sigd_yf = std(diff(yfor));
    sigd_s1 = std(diff(seg1));
    if N2>lseg
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
if N2>lseg
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
        % model2 = armax_par(seg2n,ord);
        model2 = armax(seg2n,ord);
    catch E
        msg = getReport(E);
        go = false;
        return
    end
    data2 = iddata( flipud(seg2n(2:end)), [] );
    yback = forecast( model2, data2, N2+1, opt);
    yback = yback.y;
    yback = flipud(yback);
    
    % Estimation of residuals
    e2 = resid(model2, data2);
    rstd2 = std(e2.OutputData);
    r2 = rstd2*randn(size(yback));
    yback = ( (yback + r2).*sig_s2 + mean(seg2) );
    
    sigd_yb = std(diff(yback));
    sigd_s2 = std(diff(seg2));
    if N2>lseg
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
