%% Generate Fits for Figure 5, 6 and 7
% These fits find the best fit model for a given subject group across the
% two late phases of the experiment (Late Learning, and Early Washout).

% In this test, we fix the amount learned and the the threshold for an 
% acceptable fit. This is used to develop confidence intervals for the best
% fit data. That is, finding the range of learning that is acceptable, and
% how stringent the confidence intervals can be before a solution can no
% longer be found. 

% Add supporting function paths
addpath('../../Scripts/Supporting Functions')
addpath('../Model-Specific Supporting Functions')

% For which age group
ages = {'Y','O'};
% Number of Model Fits per Proportion Learned
n = 1;
% Explore 40% - 120% Proportion Learned
alphas = 0.4:0.05:1.2;
% alphas = 0.6:0.1:0.7;

% Set confidence interval 
ci_to_test = [1.96,...  % 95% CI
              1.645,... % 90% CI
              1.281,... % 80% CI
              1.036,... % 70% CI
              0.842,... % 60% CI
              0.674];   % 50% CI

% Assign What Parameters to Vary for Model Fitting
FitIndices = [1, ... % R(1,1) = R(2,2) (control effort) (second deriv of force)
              1, ... % Q(5,5) = Q(6,6) (force penalty)
              1, ... % Q(1,1) = Q(7,7) = -Q(1,7) = -Q(7,1) (horizontal tracking)             
              1, ... % Q(2,2) = Q(8,8) = -Q(2,8) = -Q(8,2) (vertical dist. penalty)             
              0, ... % Q(3,3) % Horizontal Velocity
              1, ... % Phi(1,1) = Phi(2,2) = Phi(7,7) = Phi(8,8) = -Phi(1,7)... (terminal position)
              1, ... % Phi(3,3) = Phi(4,4) (terminal velocity)
              1, ... % Phi(5,5) = Phi(6,6) (terminal force)
              0, ... % estGain - LN1 (proportional estimate of curl field)
              0, ... % estGain - EF1 (proportional estimate of curl field)
              0, ... % estGain - LF2 (proportional estimate of curl field)
              1, ... % estGain - EN2 (proportional estimate of curl field)
              1, ... % Q(9,9) = Q(10,10) (rate of force)
              1];    % Phi(9,9) = Phi(10,10) (terminal rate of force)

% Run Fit Code (looking at learning metrics for LF2 and EN2)
for n_ci = 1:length(ci_to_test)
    ci_scale = ci_to_test(n_ci); 
    fit_model_to_group_data
    disp(['% % %-------- done with ',num2str(ci_to_test(n_ci)),'-------% % %'])
end

% Save Data