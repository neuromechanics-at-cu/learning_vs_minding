function RMS = qualityOfModelFit_AllPhases_individual(X,subjData,subjDataChan,FitIndices,Params,subjNum)


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
    Params.estGainLN1 = X(count);
    count = count + 1;
end
if FitIndices(10)
    Params.estGainEF1 = X(count);
    count = count + 1;
end
if FitIndices(11)
    Params.estGainLF2 = X(count);
    count = count + 1;
end
if FitIndices(12)
    Params.estGainEN2 = X(count);
    count = count + 1;
end
if FitIndices(13) 
    Params.Q(9,9) = X(count);
    Params.Q(10,10) = X(count);
    count = count + 1;
end
if FitIndices(14) 
    Params.Phi(9,9) = X(count);
    Params.Phi(10,10) = X(count);
end

% LN1
LN1Params = Params;
if FitIndices(9)
    LN1Params.estGain = Params.estGainLN1;
else
    LN1Params.estGain = 0;
end
LN1Params.curlGain = 0;
LN1Params.trialLength = find(~isnan(subjData.LN1.PxSubjNorm(:,subjNum)),1,'last');
forceInit_LN1 = Params.Mass*[subjData.LN1.AxSubjNorm(1,subjNum); 
                            subjData.LN1.AySubjNorm(1,subjNum)] - ...
                LN1Params.curlGain*[0 1;-1 0]*[subjData.LN1.VxSubjNorm(1,subjNum);
                                               subjData.LN1.VySubjNorm(1,subjNum)];
LN1Params.x0 = [subjData.LN1.PxSubjNorm(1,subjNum);
            subjData.LN1.PySubjNorm(1,subjNum);
            subjData.LN1.VxSubjNorm(1,subjNum);
            subjData.LN1.VySubjNorm(1,subjNum);
            forceInit_LN1(1);forceInit_LN1(2);
            0;
            0.1;0;0];
SOL_LN1 = generateReachTrajectory_CurlOrBaseline(LN1Params);
AccMat = inv(Params.Mass)*([SOL_LN1.x(:,5)';SOL_LN1.x(:,6)'] + LN1Params.curlGain*[0 1;-1 0]*[SOL_LN1.x(:,3)';SOL_LN1.x(:,4)']);
SOL_LN1.Ax = AccMat(1,:)';
SOL_LN1.Ay = AccMat(2,:)';
% LN1 CHANNEL
LN1ChanParams = Params;
if FitIndices(9)
    LN1ChanParams.estGain = Params.estGainLN1;
else
    LN1ChanParams.estGain = 0;
end
LN1ChanParams.curlGain = 0;
LN1ChanParams.trialLength = find(~isnan(subjDataChan.LN1.PxSubjNorm(:,subjNum)),1,'last');
LN1ChanParams.x0 = [subjDataChan.LN1.PxSubjNorm(1,subjNum);
                    subjDataChan.LN1.PySubjNorm(1,subjNum);
                    subjDataChan.LN1.VxSubjNorm(1,subjNum);
                    subjDataChan.LN1.VySubjNorm(1,subjNum);
                    subjDataChan.LN1.FxSubjNorm(1,subjNum);
                    subjDataChan.LN1.FySubjNorm(1,subjNum);
                    0;
                    0.1;0;0];
SOL_LN1Chan = generateReachTrajectory_Channel(LN1ChanParams);
ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_LN1Chan.x(:,1)';SOL_LN1Chan.x(:,2)';SOL_LN1Chan.x(:,3)';SOL_LN1Chan.x(:,4)'];
SOL_LN1Chan.Fx = ForceMat(1,:)';
SOL_LN1Chan.Fy = ForceMat(2,:)';
SOL_LN1Chan.AdaptInd= fitgain(SOL_LN1Chan.Fx,SOL_LN1Chan.x(:,4)); 
[~,t_ind] = max(abs(SOL_LN1Chan.Fx));
SOL_LN1Chan.FxMax = SOL_LN1Chan.Fx(t_ind);

% EF1
EF1Params = Params;
if FitIndices(10)
    EF1Params.estGain = Params.estGainEF1;
else
    EF1Params.estGain = Params.estGain;
end
EF1Params.curlGain = -20;
EF1Params.trialLength = find(~isnan(subjData.EF1.PxSubjNorm(:,subjNum)),1,'last');
forceInit_EF1 = Params.Mass*[subjData.EF1.AxSubjNorm(1,subjNum); 
                            subjData.EF1.AySubjNorm(1,subjNum)] - ...
                EF1Params.curlGain*[0 1;-1 0]*[subjData.EF1.VxSubjNorm(1,subjNum);
                                               subjData.EF1.VySubjNorm(1,subjNum)];
EF1Params.x0 = [subjData.EF1.PxSubjNorm(1,subjNum);
                subjData.EF1.PySubjNorm(1,subjNum);
                subjData.EF1.VxSubjNorm(1,subjNum);
                subjData.EF1.VySubjNorm(1,subjNum);
                forceInit_EF1(1);forceInit_EF1(2);
                0;
                0.1;0;0];
SOL_EF1 = generateReachTrajectory_CurlOrBaseline(EF1Params);
AccMat = inv(Params.Mass)*([SOL_EF1.x(:,5)';SOL_EF1.x(:,6)'] + EF1Params.curlGain*[0 1;-1 0]*[SOL_EF1.x(:,3)';SOL_EF1.x(:,4)']);
SOL_EF1.Ax = AccMat(1,:)';
SOL_EF1.Ay = AccMat(2,:)';
% EF1 CHANNEL
EF1ChanParams = Params;
if FitIndices(10)
    EF1ChanParams.estGain = Params.estGainEF1;
else
    EF1ChanParams.estGain = Params.estGain;
end
EF1ChanParams.curlGain = -20;
EF1ChanParams.trialLength = find(~isnan(subjData.EF1.PxSubjNorm(:,subjNum)),1,'last');
EF1ChanParams.x0 = [subjDataChan.EF1.PxSubjNorm(1,subjNum);
                    subjDataChan.EF1.PySubjNorm(1,subjNum);
                    subjDataChan.EF1.VxSubjNorm(1,subjNum);
                    subjDataChan.EF1.VySubjNorm(1,subjNum);
                    subjDataChan.EF1.FxSubjNorm(1,subjNum);
                    subjDataChan.EF1.FySubjNorm(1,subjNum);
                    0;
                    0.1;0;0];
SOL_EF1Chan = generateReachTrajectory_Channel(EF1ChanParams);
ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_EF1Chan.x(:,1)';SOL_EF1Chan.x(:,2)';SOL_EF1Chan.x(:,3)';SOL_EF1Chan.x(:,4)'];
SOL_EF1Chan.Fx = ForceMat(1,:)';
SOL_EF1Chan.Fy = ForceMat(2,:)';
SOL_EF1Chan.AdaptInd= fitgain(SOL_EF1Chan.Fx,SOL_EF1Chan.x(:,4)); 
[~,t_ind] = max(abs(SOL_EF1Chan.Fx));
SOL_EF1Chan.FxMax = SOL_EF1Chan.Fx(t_ind);

% LF2
LF2Params = Params;
if FitIndices(11)
    LF2Params.estGain = Params.estGainLF2;
else
    LF2Params.estGain = Params.estGain;    
end
LF2Params.curlGain = -20;
LF2Params.trialLength = find(~isnan(subjData.LF2.PxSubjNorm(:,subjNum)),1,'last');
forceInit_LF2 = Params.Mass*[subjData.LF2.AxSubjNorm(1,subjNum); 
                            subjData.LF2.AySubjNorm(1,subjNum)] - ...
                LF2Params.curlGain*[0 1;-1 0]*[subjData.LF2.VxSubjNorm(1,subjNum);
                                               subjData.LF2.VySubjNorm(1,subjNum)];
LF2Params.x0 = [subjData.LF2.PxSubjNorm(1,subjNum);
            subjData.LF2.PySubjNorm(1,subjNum);
            subjData.LF2.VxSubjNorm(1,subjNum);
            subjData.LF2.VySubjNorm(1,subjNum);
            forceInit_LF2(1);forceInit_LF2(2);
            0;
            0.1;
            0;0];
SOL_LF2 = generateReachTrajectory_CurlOrBaseline(LF2Params);
AccMat = inv(Params.Mass)*([SOL_LF2.x(:,5)';SOL_LF2.x(:,6)'] + LF2Params.curlGain*[0 1;-1 0]*[SOL_LF2.x(:,3)';SOL_LF2.x(:,4)']);
SOL_LF2.Ax = AccMat(1,:)';
SOL_LF2.Ay = AccMat(2,:)';
% LF2 CHANNEL
LF2ChanParams = Params;
if FitIndices(11)
    LF2ChanParams.estGain = Params.estGainLF2;
else
    LF2ChanParams.estGain = Params.estGain;
end
LF2ChanParams.curlGain = -20;
LF2ChanParams.trialLength = find(~isnan(subjDataChan.LF2.PxSubjNorm(:,subjNum)),1,'last');
LF2ChanParams.x0 = [subjDataChan.LF2.PxSubjNorm(1,subjNum);
                    subjDataChan.LF2.PySubjNorm(1,subjNum);
                    subjDataChan.LF2.VxSubjNorm(1,subjNum);
                    subjDataChan.LF2.VySubjNorm(1,subjNum);
                    subjDataChan.LF2.FxSubjNorm(1,subjNum);
                    subjDataChan.LF2.FySubjNorm(1,subjNum);
                    0;
                    0.1;0;0];
SOL_LF2Chan = generateReachTrajectory_Channel(LF2ChanParams);
ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_LF2Chan.x(:,1)';SOL_LF2Chan.x(:,2)';SOL_LF2Chan.x(:,3)';SOL_LF2Chan.x(:,4)'];
SOL_LF2Chan.Fx = ForceMat(1,:)';
SOL_LF2Chan.Fy = ForceMat(2,:)';
SOL_LF2Chan.AdaptInd= fitgain(SOL_LF2Chan.Fx,SOL_LF2Chan.x(:,4)); 
[~,t_ind] = max(abs(SOL_LF2Chan.Fx));
SOL_LF2Chan.FxMax = SOL_LF2Chan.Fx(t_ind);

% EN2
EN2Params = Params;
if FitIndices(12)
    EN2Params.estGain = Params.estGainEN2;
else
    EN2Params.estGain = subjData.EN2.estGainMean;
end
EN2Params.curlGain = 0;
EN2Params.trialLength = find(~isnan(subjData.EN2.PxSubjNorm(:,subjNum)),1,'last');
forceInit_EN2 = Params.Mass*[subjData.EN2.AxSubjNorm(1,subjNum); 
                            subjData.EN2.AySubjNorm(1,subjNum)] - ...
                EN2Params.curlGain*[0 1;-1 0]*[subjData.EN2.VxSubjNorm(1,subjNum);
                                               subjData.EN2.VySubjNorm(1,subjNum)];
EN2Params.x0 = [subjData.EN2.PxSubjNorm(1,subjNum);
            subjData.EN2.PySubjNorm(1,subjNum);
            subjData.EN2.VxSubjNorm(1,subjNum);
            subjData.EN2.VySubjNorm(1,subjNum);
            forceInit_EN2(1);forceInit_EN2(2);
            0;
            0.1;
            0;0];
SOL_EN2 = generateReachTrajectory_CurlOrBaseline(EN2Params);
AccMat = inv(Params.Mass)*([SOL_EN2.x(:,5)';SOL_EN2.x(:,6)'] + EN2Params.curlGain*[0 1;-1 0]*[SOL_EN2.x(:,3)';SOL_EN2.x(:,4)']);
SOL_EN2.Ax = AccMat(1,:)';
SOL_EN2.Ay = AccMat(2,:)';
% EN2 CHANNEL
EN2ChanParams = Params;
if FitIndices(12)
    EN2ChanParams.estGain = Params.estGainEN2;
else
    EN2ChanParams.estGain = subjDataChan.EN2.estGainMean;
end
EN2ChanParams.curlGain = 0;
EN2ChanParams.trialLength = find(~isnan(subjDataChan.EN2.PxSubjNorm(:,subjNum)),1,'last');
EN2ChanParams.x0 = [subjDataChan.EN2.PxSubjNorm(1,subjNum);
                    subjDataChan.EN2.PySubjNorm(1,subjNum);
                    subjDataChan.EN2.VxSubjNorm(1,subjNum);
                    subjDataChan.EN2.VySubjNorm(1,subjNum);
                    subjDataChan.EN2.FxSubjNorm(1,subjNum);
                    subjDataChan.EN2.FySubjNorm(1,subjNum);
                    0;
                    0.1;0;0];
SOL_EN2Chan = generateReachTrajectory_Channel(EN2ChanParams);
ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_EN2Chan.x(:,1)';SOL_EN2Chan.x(:,2)';SOL_EN2Chan.x(:,3)';SOL_EN2Chan.x(:,4)'];
SOL_EN2Chan.Fx = ForceMat(1,:)';
SOL_EN2Chan.Fy = ForceMat(2,:)';
SOL_EN2Chan.AdaptInd= fitgain(SOL_EN2Chan.Fx,SOL_EN2Chan.x(:,4)); 
[~,t_ind] = max(abs(SOL_EN2Chan.Fx));
SOL_EN2Chan.FxMax = SOL_EN2Chan.Fx(t_ind);

%% Calculate Objective Function

%LN1
[~,mpx_ind] = max(abs(subjData.LN1.PxSubjNorm(:,subjNum)));
[~,mpx_ind_mdl] = max(abs(SOL_LN1.x(:,1)));
endpt = find(~isnan(subjData.LN1.PxSubjNorm(:,subjNum)),1,'last');

if isnan(subjDataChan.LN1.FxMax(subjNum))
    RMS_LN1 = (subjData.LN1.PxSubjNorm(endpt,subjNum) - SOL_LN1.x(end,1))^2 + ... % End Px
          (subjData.LN1.PySubjNorm(endpt,subjNum) - SOL_LN1.x(end,2))^2 + ... % End Py
          10*(subjData.LN1.PxSubjNorm(mpx_ind,subjNum) - SOL_LN1.x(mpx_ind_mdl,1))^2 + ... % Max Px
          10*(subjData.LN1.PySubjNorm(mpx_ind,subjNum) - SOL_LN1.x(mpx_ind_mdl,2))^2;

else
    RMS_LN1 = (subjData.LN1.PxSubjNorm(endpt,subjNum) - SOL_LN1.x(end,1))^2 + ... % End Px
          (subjData.LN1.PySubjNorm(endpt,subjNum) - SOL_LN1.x(end,2))^2 + ... % End Py
          10*(subjData.LN1.PxSubjNorm(mpx_ind,subjNum) - SOL_LN1.x(mpx_ind_mdl,1))^2 + ... % Max Px
          10*(subjData.LN1.PySubjNorm(mpx_ind,subjNum) - SOL_LN1.x(mpx_ind_mdl,2))^2 + ... % Py of Max Px
          (subjDataChan.LN1.FxMax(subjNum) - SOL_LN1Chan.FxMax)^2 + ... % Max Fx
          (subjDataChan.LN1.estGain(subjNum) - SOL_LN1Chan.AdaptInd)^2; % Adapt Ind
end

%EF1
[~,mpx_ind] = max(abs(subjData.EF1.PxSubjNorm(:,subjNum)));
[~,mpx_ind_mdl] = max(abs(SOL_EF1.x(:,1)));
endpt = find(~isnan(subjData.EF1.PxSubjNorm(:,subjNum)),1,'last');

if isnan(subjDataChan.EF1.FxMax(subjNum))
    RMS_EF1 = (subjData.EF1.PxSubjNorm(endpt,subjNum) - SOL_EF1.x(end,1))^2 + ... % End Px
          (subjData.EF1.PySubjNorm(endpt,subjNum) - SOL_EF1.x(end,2))^2 + ... % End Py
          10*(subjData.EF1.PxSubjNorm(mpx_ind,subjNum) - SOL_EF1.x(mpx_ind_mdl,1))^2 + ... % Max Px
          10*(subjData.EF1.PySubjNorm(mpx_ind,subjNum) - SOL_EF1.x(mpx_ind_mdl,2))^2;

else
    RMS_EF1 = (subjData.EF1.PxSubjNorm(endpt,subjNum) - SOL_EF1.x(end,1))^2 + ... % End Px
          (subjData.EF1.PySubjNorm(endpt,subjNum) - SOL_EF1.x(end,2))^2 + ... % End Py
          10*(subjData.EF1.PxSubjNorm(mpx_ind,subjNum) - SOL_EF1.x(mpx_ind_mdl,1))^2 + ... % Max Px
          10*(subjData.EF1.PySubjNorm(mpx_ind,subjNum) - SOL_EF1.x(mpx_ind_mdl,2))^2 + ... % Py of Max Px
          (subjDataChan.EF1.FxMax(subjNum) - SOL_EF1Chan.FxMax)^2 + ... % Max Fx
          (subjDataChan.EF1.estGain(subjNum) - SOL_EF1Chan.AdaptInd)^2; % Adapt Ind
end

%LF2
[~,mpx_ind] = max(abs(subjData.LF2.PxSubjNorm(:,subjNum)));
[~,mpx_ind_mdl] = max(abs(SOL_LF2.x(:,1)));
endpt = find(~isnan(subjData.LF2.PxSubjNorm(:,subjNum)),1,'last');

if isnan(subjDataChan.LF2.FxMax(subjNum))
    RMS_LF2 = (subjData.LF2.PxSubjNorm(endpt,subjNum) - SOL_LF2.x(end,1))^2 + ... % End Px
          (subjData.LF2.PySubjNorm(endpt,subjNum) - SOL_LF2.x(end,2))^2 + ... % End Py
          10*(subjData.LF2.PxSubjNorm(mpx_ind,subjNum) - SOL_LF2.x(mpx_ind_mdl,1))^2 + ... % Max Px
          10*(subjData.LF2.PySubjNorm(mpx_ind,subjNum) - SOL_LF2.x(mpx_ind_mdl,2))^2;

else
    RMS_LF2 = (subjData.LF2.PxSubjNorm(endpt,subjNum) - SOL_LF2.x(end,1))^2 + ... % End Px
          (subjData.LF2.PySubjNorm(endpt,subjNum) - SOL_LF2.x(end,2))^2 + ... % End Py
          10*(subjData.LF2.PxSubjNorm(mpx_ind,subjNum) - SOL_LF2.x(mpx_ind_mdl,1))^2 + ... % Max Px
          10*(subjData.LF2.PySubjNorm(mpx_ind,subjNum) - SOL_LF2.x(mpx_ind_mdl,2))^2 + ... % Py of Max Px
          (subjDataChan.LF2.FxMax(subjNum) - SOL_LF2Chan.FxMax)^2 + ... % Max Fx
          (subjDataChan.LF2.estGain(subjNum) - SOL_LF2Chan.AdaptInd)^2; % Adapt Ind
end

%EN2
[~,mpx_ind] = max(abs(subjData.EN2.PxSubjNorm(:,subjNum)));
[~,mpx_ind_mdl] = max(abs(SOL_EN2.x(:,1)));
endpt = find(~isnan(subjData.EN2.PxSubjNorm(:,subjNum)),1,'last');

if isnan(subjDataChan.EN2.FxMax(subjNum))
    RMS_EN2 = (subjData.EN2.PxSubjNorm(endpt,subjNum) - SOL_EN2.x(end,1))^2 + ... % End Px
          (subjData.EN2.PySubjNorm(endpt,subjNum) - SOL_EN2.x(end,2))^2 + ... % End Py
          10*(subjData.EN2.PxSubjNorm(mpx_ind,subjNum) - SOL_EN2.x(mpx_ind_mdl,1))^2 + ... % Max Px
          10*(subjData.EN2.PySubjNorm(mpx_ind,subjNum) - SOL_EN2.x(mpx_ind_mdl,2))^2;

else
    RMS_EN2 = (subjData.EN2.PxSubjNorm(endpt,subjNum) - SOL_EN2.x(end,1))^2 + ... % End Px
          (subjData.EN2.PySubjNorm(endpt,subjNum) - SOL_EN2.x(end,2))^2 + ... % End Py
          10*(subjData.EN2.PxSubjNorm(mpx_ind,subjNum) - SOL_EN2.x(mpx_ind_mdl,1))^2 + ... % Max Px
          10*(subjData.EN2.PySubjNorm(mpx_ind,subjNum) - SOL_EN2.x(mpx_ind_mdl,2))^2 + ... % Py of Max Px
          (subjDataChan.EN2.FxMax(subjNum) - SOL_EN2Chan.FxMax)^2 + ... % Max Fx
          (subjDataChan.EN2.estGain(subjNum) - SOL_EN2Chan.AdaptInd)^2; % Adapt Ind
end


% Combine Phases  
RMS = RMS_LN1 + RMS_EF1 + RMS_LF2 + RMS_EN2;