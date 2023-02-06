%% Generate Model on Model Fits - Fig 3

% add supporting fuctions path
addpath('../../Scripts/Supporting Functions')
addpath('../Model-Specific Supporting Functions')

%% Define Model-Generated Trajectory To Be Fit
% Cost Ratio Type
cost_ratio_type = 1; %1 - high, 2 - med, 3 - low

% Proportion Learned
prop_learned = 1.0;
% prop_learned = 0.6;

% Which Trajectories
traj_type = 1;% 1 - learning only, 2- learning and washout

%% Define Search Parameters
% What parameters will be allowed to vary 
vary_learning = 1; %1 = true, 0 = false
vary_costs = 1;  %1 = true, 0 = false

% How many fits
n_sols = 10;
% Starting alpha
alphas = 0.7;
% Max Objective Function
max_objective_fxn = 1e-3;

% If learning isnt varied, search start needs to match true learning
if vary_learning == 0
    alphas = prop_learned;
end

%% Run Fit Code (looking at learning metrics for LF2 and EN2)
if traj_type == 1
    fit_model_to_model_learningonly
elseif traj_type == 2
    fit_model_to_model_learningandwashout
end
