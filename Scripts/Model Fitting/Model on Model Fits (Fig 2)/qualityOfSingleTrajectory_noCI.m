%% Calculate Quality of Model Fit
% This function is what is minimized 
% (smaller number ~ better fitting model)
%
%
% Data from both Channel Trials and Regular trials are used
% 
% Quality of Fit looks at how closely the position of maximum perpendicular
% error matches, the end point, the maximum perpendicular force, and
% adaptation index. 


function [RMS] = qualityOfSingleTrajectory_noCI(X,SOL,SOL_Chan,FitIndices,Params)

stiffConst = 2000;
dampConst = 50;

count = 1;
if FitIndices(1)
    Params.R(1,1) = X(count);
    Params.R(2,2) = X(count);
    count = count + 1;
end
if FitIndices(2)
    Params.Q(5,5) = X(count);
    Params.Q(6,6) = X(count);
    count = count + 1;
end
if FitIndices(3)
    Params.Q(1,1) = X(count);
    Params.Q(7,7) = X(count);
    Params.Q(1,7) = -1*X(count);
    Params.Q(7,1) = -1*X(count);
    count = count + 1;
end
if FitIndices(4) 
    Params.Q(2,2) = X(count);
    Params.Q(8,8) = X(count);
    Params.Q(2,8) = -1*X(count);
    Params.Q(8,2) = -1*X(count);
    count = count + 1;
end
if FitIndices(5)
    Params.Q(3,3) = X(count);
    count = count + 1;
end
if FitIndices(6) 
    Params.Phi(1,1) = X(count);
    Params.Phi(7,7) = X(count);
    Params.Phi(1,7) = -1*X(count);
    Params.Phi(7,1) = -1*X(count);
    Params.Phi(2,2) = X(count);
    Params.Phi(8,8) = X(count);
    Params.Phi(2,8) = -1*X(count);
    Params.Phi(8,2) = -1*X(count);
    count = count + 1;
end
if FitIndices(7)
    Params.Phi(3,3) = X(count);
    Params.Phi(4,4) = X(count);
    count = count + 1;
end
if FitIndices(8) 
    Params.Phi(5,5) = X(count);
    Params.Phi(6,6) = X(count);
    count = count + 1;
end
if FitIndices(9)
    Params.estGain = X(count);
    count = count + 1;
end
if FitIndices(10) 
    Params.Q(9,9) = X(count);
    Params.Q(10,10) = X(count);
    count = count + 1;
end
if FitIndices(11) 
    Params.Phi(9,9) = X(count);
    Params.Phi(10,10) = X(count);
end

%% Generate Model-Produced Trajectories

Params.trialLength = length(SOL.x(:,1));
% Params.x0 = 

SOL_new = generateReachTrajectory_CurlOrBaseline(Params);

AccMat = inv(Params.Mass)*([SOL_new.x(:,5)';SOL_new.x(:,6)'] + Params.curlGain*[0 1;-1 0]*[SOL_new.x(:,3)';SOL_new.x(:,4)']);
SOL_new.Ax = AccMat(1,:)';
SOL_new.Ay = AccMat(2,:)';
% CHANNEL
ChanParams = Params;
ChanParams.curlGain = -20;
ChanParams.trialLength = length(SOL_Chan.x(:,1));
ChanParams.x0 = SOL_Chan.x(1,:);
SOL_Chan_new = generateReachTrajectory_Channel(ChanParams);

ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_Chan_new.x(:,1)';SOL_Chan_new.x(:,2)';SOL_Chan_new.x(:,3)';SOL_Chan_new.x(:,4)'];
SOL_Chan_new.Fx = ForceMat(1,:)';
SOL_Chan_new.Fy = ForceMat(2,:)';
SOL_Chan_new.AdaptInd= fitgain(SOL_Chan_new.Fx,SOL_Chan_new.x(:,4)); 
[~,t_ind] = max(abs(SOL_Chan_new.Fx));
SOL_Chan_new.FxMax = SOL_Chan_new.Fx(t_ind);


%% Calculate Objective Function

%LF2
[~,mpx_ind] = max(abs(SOL.x(:,1)));
[~,mpx_ind_new] = max(abs(SOL_new.x(:,1)));


% RMS = (SOL_new.x(end,1) - SOL.x(end,1))^2 + ... % End Px
%       (SOL_new.x(end,2) - SOL.x(end,2))^2 + ... % End Py
%       10*(SOL_new.x(mpx_ind_new,1) - SOL.x(mpx_ind,1))^2 + ... % Max Px
%       10*(SOL_new.x(mpx_ind_new,2) - SOL.x(mpx_ind,2))^2 + ... % Py of Max Px
%       (SOL_Chan_new.FxMax - SOL_Chan.FxMax)^2 + ... % Max Fx
%       (SOL_Chan_new.AdaptInd - SOL_Chan.AdaptInd)^2; % Adapt Ind

% Experiment - 20AUG2022
RMS = (SOL_new.x(end,1) - SOL.x(end,1))^2 + ... % End Px
      (SOL_new.x(end,2) - SOL.x(end,2))^2 + ... % End Py
      2*(SOL_new.x(mpx_ind_new,1) - SOL.x(mpx_ind,1))^2 + ... % Max Px
      2*(SOL_new.x(mpx_ind_new,2) - SOL.x(mpx_ind,2))^2 + ... % Py of Max Px
      0.1*(SOL_Chan_new.FxMax - SOL_Chan.FxMax)^2 + ... % Max Fx
      0.1*(SOL_Chan_new.AdaptInd - SOL_Chan.AdaptInd)^2; % Adapt Ind

% RMS = (SOL_new.x(end,1) - SOL.x(end,1))^2 + ... % End Px
%       (SOL_new.x(end,2) - SOL.x(end,2))^2 + ... % End Py
%       10*(SOL_new.x(mpx_ind_new,1) - SOL.x(mpx_ind,1))^2 + ... % Max Px
%       10*(SOL_new.x(mpx_ind_new,2) - SOL.x(mpx_ind,2))^2; % Py of Max Px

