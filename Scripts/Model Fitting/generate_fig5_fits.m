%% Generate Fits for Figure 5
%
ages = {'Y','O'};
% Ten Model Fits per Proportion Learned
n = 10;
% Explore Intersecting Proportion Learned Between Old and Young
alphas = 0.6:0.05:0.8;
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
fit_model_to_data

% Process Cost Weight Data
costratios_from_fits
