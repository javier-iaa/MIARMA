function ic = aicplus(model, criterion)
%AICPLUS Computes Information Criteria(AIC) from a model
%   IC = AICPLUS(Model, criterion)
%
%   Model = Any IDMODEL or IDNLMODEL
%
%   AIC = Akaikes Information Criterion log(V) + 2d/N
%   where V is the loss function, d is the number of estimated parameters
%   and N is the number of estimation data.
%
%   criterion - model selection criteria, choose between:
%      'AIC'(default): Akaike Information Criterion
%      'AICc': AIC with correction for small samples
%      'BIC': Bayesian Information Criterion
%      'FPE': Final Prediction Error
%      'HQ': Hannan-Quinn Criterion
%
%   Author: Javier Pascual-Grandao
%   Revision: 1.0  
%   Date: 2019/12/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    criterion = 'AIC'; % Set default criterion
end

d = sum(model.na) + sum(model.nc) + sum(model.nb); % number of parameters
k = d + 1;
es = pvget(model,'EstimationInfo');
V = es.LossFcn;
FPE = es.FPE;
N = es.DataLength;
Nnpar = FPE/V - 1;  % = 2d/N when d<<N, i.e. the Akaike asymptotic approx.

switch criterion
    % Akaike Information Criterion
    case 'AIC'
    	ic = log(V) + Nnpar; % -2*log(L) + 2*(p+q+1)
    
    % AIC with correction for small samples
    case 'AICc'
        ic = log(V) + Nnpar + (1/N)*2*k*(k+1)/(N-k-1);
    
    % Bayesian Information Criterion (also known as the Schwarz-Bayesfunction ic = aicplus(model, criterion)
    % Criterion) -> to be revised
    case 'BIC'
        ic = log(V/N) + Nnpar*log(N)/2;

    % Akaike Final Prediction Error
    case 'FPE'    
        ic = FPE;
        
    % Hannan-Quinn Criterion -> to be revised
    case 'HQ' 
        ic = log(V) + Nnpar*log(log(N));
        
    otherwise
        error('Selection criterion is not defined properly')
end

end
