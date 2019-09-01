function YP = pred(model,data,K, Init)
% FORECAST Forecast a time series K steps into the future
%
% YP = PRED(MODEL,DATA,K)
%
% DATA: Existing data up to time N, an IDDATA object.
% MODEL: The model as any IDMODEL object, IDPOLY, IDSS, IDARX or IDGREY.
% K: The time horizon of forecasting, a positive integer with the number of samples
% YP: The forecasted output after time N, an IDDATA object with output
% only, covering the time span N+1:N+K.
%
%  YP = FORECAST(MODEL,DATA,K, INIT)
%  wehere INIT is 'z' or 'e' allows specification of initial conditions (at
%  time = Data.SamplingInstants(1)).
%
% See also idmodel/predict, which computes a fixed horizon prediction
% along the existing data record.
%
% Based on the function posted online by Rajiv Singh on 22 Jul 2011 in
% matlabcentral
% By Javier Pascual-Granado, 2014
% <a href="matlab:web http://www.iaa.es;">IAA-CSIC, Spain</a>

% Version: 1.0.1
% Changes from last version: x0est function is deprecated, now we use
%   findstates to estimate the initial state when init == 'e'
% Last update: 1-sep-2015

[N, ny] = size(data); % assume data is iddata

Mss = idss(model);

ord = size(pvget(Mss,'A'),1);

if ord>N
 error('Forecast:TooFewSamples','The data should contain at least %d samples.',ord)
end

% Default
if nargin<4, 
    Init = 'e'; 
end

yp = zeros(K,ny);
mp = getPredictor(Mss);
[Ap,Bp,Cp] = ssdata(mp);

if Init=='z',
    xt = ltitr(Ap, Bp, data.y); % use zero init
    x0 = xt(end,:);
else % Init == 'e'
%     [A1,B1,C1,D1,K1] = ssdata(Mss);
    x00 = findstates(model,data);
%     x00 = x0est(data.y,A1,B1,C1,D1,K1,size(C1,1),size(B1,2),...
%         250e3,eye(size(C1,1)));
    x0 = ltitr(Ap,Bp,data.y,x00); x0 = x0(end,:);
end

u = [data.y(end,:); zeros(1,ny)];
for ct = 1:K
    xt = ltitr(Ap, Bp, u, x0);
    x0 = xt(end,:);
    yp(ct,:) = (Cp*x0.').';
    u = [yp(ct,:); zeros(1,ny)];
end

YP = data; YP.y = yp; YP.u = []; YP.Name = '';
YP.UserData = []; YP.Notes = '';
YP.Tstart = data.Ts*(N+1);

%------local function -----------------------------------
function mp = getPredictor(sysd)

[A,B,C,D,K] = ssdata(sysd);
Ts = sysd.Ts;
Ny = size(D,1);
mp = idss(A-K*C, [K B-K*D], C, [zeros(Ny), D], zeros(size(A,1),Ny), 'Ts', Ts);
mp.InputDelay = [zeros(Ny,1); sysd.InputDelay];
