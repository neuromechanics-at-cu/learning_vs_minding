%% Generate Fits for Figure 3
%
ages = {'Y','O'};
% One Model Fit per Age Group and Proportion Learned
n = 1;
% Set initial alpha value of 0.8 (this is arbitrary and will change
alphas = 0.8; 

% Assign What Parameters to Vary for Model Fitting
FitIndices = [1, ... % R(1,1) = R(2,2) (control effort) (second deriv of force)
              1, ... % Q(5,5) = Q(6,6) (force penalty)
              1, ... % Q(1,1) = Q(7,7) = -Q(1,7) = -Q(7,1) (horizontal tracking)             
              1, ... % Q(2,2) = Q(8,8) = -Q(2,8) = -Q(8,2) (vertical dist. penalty)             
              0, ... % Q(3,3) % Horizontal Velocity
              1, ... % Phi(1,1) = Phi(2,2) = Phi(7,7) = Phi(8,8) = -Phi(1,7)... (terminal position)
              1, ... % Phi(3,3) = Phi(4,4) (terminal velocity)
              1, ... % Phi(5,5) = Phi(6,6) (terminal force)
              1, ... % estGain - LN1 (proportional estimate of curl field)
              1, ... % estGain - EF1 (proportional estimate of curl field)
              1, ... % estGain - LF2 (proportional estimate of curl field)
              1, ... % estGain - EN2 (proportional estimate of curl field)
              1, ... % Q(9,9) = Q(10,10) (rate of force)
              1];    % Phi(9,9) = Phi(10,10) (terminal rate of force)

% Run Fit Code
fit_model_to_data
