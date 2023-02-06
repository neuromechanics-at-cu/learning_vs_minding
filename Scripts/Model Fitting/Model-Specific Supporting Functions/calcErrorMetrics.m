function metrics = calcErrorMetrics(Params)

% Generate Trajectory
SOL = generateReachTrajectory_CurlOrBaseline(Params);
SOL_Chan = generateReachTrajectory_Channel(Params);

% Max Px
[~,mpx_ind] = max(abs(SOL.x(:,1)));
metrics.MaxPx = SOL.x(mpx_ind,1);

% Max Fx
stiffConst = 2000;
dampConst = 50;
ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_Chan.x(:,1)';SOL_Chan.x(:,2)';SOL_Chan.x(:,3)';SOL_Chan.x(:,4)'];
SOL_Chan.Fx = ForceMat(1,:)';
SOL_Chan.Fy = ForceMat(2,:)';
[~,mfx_ind] = max(abs(SOL_Chan.Fx));
metrics.MaxFx = SOL_Chan.Fx(mfx_ind);

% Adapt Ind
metrics.adaptInd = fitgain(SOL_Chan.Fx,SOL_Chan.x(:,4)); 