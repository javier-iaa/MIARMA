function yi = ftcorr(y, stat, varargin)
% Fourier correction to the ARMA interpolation based on Fourier 
% reconstruction.
%
% The time series is assumed to be evenly sampled so no time vector 
% is necessary as input.
%
% Inputs:       y - observational data
%                   stat - flag to identify the gaps
%                   varargin can be: 'win', 'conf', 'ofac', 'range' followed by the 
%                       corresponding value.
%                       win - true to use a window, false otherwise (default)
%                       ofac - oversampling factor for the LS computation (default = 1)
%                       range - factor that controls the range of frequencies 
%                           to explore, e.g. r=1 for Nyquist (default)
%
% Outputs:    yf - reconstructed data
%
% Date: 24-mar-2021
% Javier Pascual-Granado

%% Parameters

% Frequency range to explore
range_flag = find( strcmp(varargin,'range') );
if ~isempty(range_flag)
    r = varargin{range_flag + 1};
else
    r = 1; % Nyquist frequency
end

% Windowing sometimes improve the estimation
win_flag = find( strcmp(varargin,'win') );
if ~isempty(win_flag)
    win = varargin{win_flag + 1};
else
    win = false;
end

% Oversampling factor for the calculation of the FFT
ofac_flag = find( strcmp(varargin,'ofac') );
if ~isempty(ofac_flag)
    ofac = varargin{ofac_flag + 1};
else
    ofac = 1;
end

% Confidence level
% conf_flag = find( strcmp(varargin,'conf') );
% if ~isempty(conf_flag)
%     conf = varargin{conf_flag + 1};
% else
%     conf = 99;
% end

% alpha = (100 - conf)/2/100;
% a = norminv( [alpha, 1-alpha], 0, 1 );

% Cutoff level
% level=a(2) ;
level = 10;

%% Prepare data
sy = size(y);

if sy(2)==1
    y = y';
end

n = length(y);

my = mean(y);

% Reconstruction is sometimes a bit better with windowing
if win
    y = y - my;
    window = nuttallwin(n);
	y = y.*window';
    y = y + my;
end

%% Calculate the periodogram and estimate noise level

nf = ofac*n;

% Normalization
% fac = sqrt((2*nf+1)/2);
% fac = (2*sqrt(2)/n);
fac = 1;

% Calculate the DFT
F = fac*fft(y, nf);
aF = abs(F);
ph = angle(F);

% Frequencies corresponding to the DFT calculation
% T = x(end) - x(1);
% df = 1/T/ofac;
% fnyq = 1/(x(2)-x(1))/2;
% wk1 = 0:df:2*fnyq;

% Consider only a part of the spectrum for confidence estimation
spec = F( 2:fix(nf/2/r) );

%  The noise level is level*mean power of spectra
noise_est = sqrt( median( spec.*conj(spec) ) );
noise_limit = sqrt( level ) * noise_est;

%% Signal extraction and clean of the amplitude spectrum

% Find the highest amplitude peak
[p, l] = findpeaks( aF, 'MINPEAKHEIGHT', noise_limit, 'MINPEAKDISTANCE', 3*ofac);

if ~isempty(p)
    F0 = F;
    % Clean the amplitude spectrum of other components than signal
    F(:) = 0.0;
    F(1) = F0(1);
    F(l) = p.*exp( 1i*ph(l) );

end

%% Inverse Fourier transform    
Fb = ifft(F)/fac;
yf = real(Fb);
% yi = imag(Fb);

% window correction (rescaling to the state without windowing)
if win
    yf = yf./window';
end

yi = y;
nz = length( find(stat==1) );
noi = noise_est*sqrt(pi/nf)*randn( 1, nz );
yi(stat==1) = yf(stat==1) + noi;

end
