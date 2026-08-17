function myalg = algoprop()
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
% myalg.InitialState = 'Auto';
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