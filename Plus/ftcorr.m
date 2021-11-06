function yi = ftcorr(y, stat, varargin)
% Fourier correction to the ARMA interpolation based on Fourier 
% reconstruction.
%
% The time series is assumed to be evenly sampled so no time vector 
% is necessary as input.
%
% Inputs:       y - observational data
%                   stat - flag to identify the gaps
%                   varargin can be: 'conf', 'ofac', 'range' followed by the 
%                       corresponding value.
%                       ofac - oversampling factor for the LS computation (default = 1)
%                       range - factor that controls the range of frequencies 
%                           to explore, e.g. r=1 for Nyquist (default)
%
% Outputs:    yf - reconstructed data
%
% Version: 0.2
% Changes:
% - Minor fixes
%
% Date: 29/10/2021
% Javier Pascual-Granado

% Prepare data
sy = size(y);

if sy(2)==1
    y = y';
end

n = length(y);

%% Parameters

% Frequency range to explore
range_flag = find( strcmp(varargin,'range') );
if ~isempty(range_flag)
    r = varargin{range_flag + 1};
else
    r = 1; % Nyquist frequency
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
level = 100;

%% Calculate the periodogram and estimate noise level

nf = ofac*n;

% Normalization
% fac = sqrt((2*nf+1)/2);
% fac = (2*sqrt(2)/n);
fac = 1;

% Calculate the DFT
F0 = fac*fft(y, nf);
aF0 = abs(F0);

% Frequencies corresponding to the DFT calculation
% T = x(end) - x(1);
% df = 1/T/ofac;
% fnyq = 1/(x(2)-x(1))/2;
% wk1 = 0:df:2*fnyq;

% Consider only a part of the spectrum for confidence estimation
spec0 = F0( 2:fix(nf/2/r) );

% Noise limit
noise_est0 = sqrt( median( spec0.*conj(spec0) ) );
noise_limit0 = sqrt( level ) * noise_est0;

%% Signal extraction and clean of the amplitude spectrum

% Find all amplitude values above the noise limit
l = find( aF0>noise_limit0 );
% [~, l] = findpeaks( aF0, 'MINPEAKHEIGHT', noise_limit0);

% Produce a new amplitude spectrum with the frequencies found
if ~isempty(l)
%     F(:) = 0.0;
%     F(:) = noise_est0;
    F = noise_est0*randn(1,nf);
    F(1) = F0(1);
    F(l) = F0(l);
else
    yi = y;
    fprintf(2, '\n  ft_corr didn''t find any significant component\n\n');
    return
end

%% Inverse Fourier transform    
Fb = ifft(F)/fac;
yf = real(Fb);

yi = y;
% nz = length( find(stat==1) );
% noi = noise_est0*sqrt(pi/nf)*randn( 1, nz );
% yi(stat==1) = yf(stat==1) + noi;
yi(stat==1) = yf(stat==1);

end
