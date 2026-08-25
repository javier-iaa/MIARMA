function inStruct = defpars(subStruct)
% function params = default_params()
% Populate a parameter structure with default values to be used in MIARMA.
% 
% Input:       subStruct - type of structure: 'params' or 'flags' for
%               parameter or flag substructure, respectively.
%
% Output:      inStruct - parameters structure
%
% Version: 0.1
%
% Author(s): Javier Pascual-Granado
% Date: 24/08/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(subStruct,'params')

    % This is the physical memory available. In case, it is different change
    % this number. Matlab R2024 is not prepared to determine this in Mac so I
    % prefer to make the program system agnostic by setting this number myself.
    inStruct.mem = 16;
    
    % Maximum length of the segment used to calculate ARMA order
    % If the optimal model does not pass the tests this will be increased until the 
    % maximum possible length
    inStruct.mseg = 1000;
    
    % Max. ratio between segment length and number of parameters for the model
    inStruct.facmax = 6;
    
    % Min. ratio between segment length and number of parameters for the model
    inStruct.facmin = 4;
    
    % Min. ratio between interpolated datapoints and the length of the segments
    inStruct.facint = 3;
    
    % Lower limit in data segment length for the ARMA interpolation.
    inStruct.npz = 36;
    
    % Lower limit in gap length in order to use ARMA interpolation, below this 
    % limit a simpler interpolation is used
    inStruct.npi = 4;
    
    % Range for the search of the optimal ARMA orders [pmin,pmax] 
    inStruct.pmin = 2;
    inStruct.pmax = 30;
    
    % The MA order is search in the range [0,qmax]
    inStruct.qmax = 30;
    
    % Parameter that set the cutoff level to extract significant frequencies with the FFT
    % that are used for FT correction in ft_corr subroutine
    inStruct.cutoff_level = 100;

elseif strcmp(subStruct,'flags')

    % Decides whether to save the Akaike matrix at a temp file
    inStruct.temp = false;
    
    % Always interpolate or not
    inStruct.always_int = true;
    
    % Flag to activate the FT correction of the arma interpolation
    inStruct.ft_corr = false;
    
    % Flag to activate (or not) ascii output with parameters and other info
    inStruct.ascii_struct = false;
    
    % Flag to activate recovery of excluded data points at the end of the 
    % gap-filling process. False means that some data points will be substituted 
    % by interpolated data points in the resulting array.
    inStruct.reco_flag = false;
    
    % Parameter flag that is used to activate the automatic search of the
    % optimal order in armaord by using an incremental Akaike matrix.
    % auto_flag is deactivated if any of pmin, pmax, qmax are also given as
    % input. If true autoarmaord.m is used instead of armaord.m
    inStruct.auto_flag = true;

    % Flag to activate debug mode and save csv files containing info about
    % the performance in gap-filling
    inStruct.debug_flag = false;

else
    error('Wrong structure or parameter.')
end