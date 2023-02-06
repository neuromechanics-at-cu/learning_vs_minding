%% Calculate Quality of Model Fit
% This function is what is minimized 
% (smaller number ~ better fitting model)
%
% Data from both Channel Trials and Regular trials are used
% 
% Quality of Fit looks at how closely the position of maximum perpendicular
% error matches, the end point, the maximum perpendicular force, and
% adaptation index. 
%
% LN1, EF1, EN2, and LF2 are equally weighted.

function [RMS] = qualityOfModelFit_AllPhases(X,subjData,subjDataChan,FitIndices,Params,ci_scale)

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
LN1Params.trialLength = length(subjData.LN1.PxSubjMean);
forceInit_LN1 = Params.Mass*[subjData.LN1.AxSubjMean(1); 
                            subjData.LN1.AySubjMean(1)] - ...
                LN1Params.curlGain*[0 1;-1 0]*[subjData.LN1.VxSubjMean(1);
                                               subjData.LN1.VySubjMean(1)];
LN1Params.x0 = [subjData.LN1.PxSubjMean(1);
            subjData.LN1.PySubjMean(1);
            subjData.LN1.VxSubjMean(1);
            subjData.LN1.VySubjMean(1);
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
LN1ChanParams.trialLength = length(subjDataChan.LN1.PxSubjMean);
LN1ChanParams.x0 = [subjDataChan.LN1.PxSubjMean(1);
                    subjDataChan.LN1.PySubjMean(1);
                    subjDataChan.LN1.VxSubjMean(1);
                    subjDataChan.LN1.VySubjMean(1);
                    subjDataChan.LN1.FxSubjMean(1);
                    subjDataChan.LN1.FySubjMean(1);
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
EF1Params.trialLength = length(subjData.EF1.PxSubjMean);
forceInit_EF1 = Params.Mass*[subjData.EF1.AxSubjMean(1); 
                            subjData.EF1.AySubjMean(1)] - ...
                EF1Params.curlGain*[0 1;-1 0]*[subjData.EF1.VxSubjMean(1);
                                               subjData.EF1.VySubjMean(1)];
EF1Params.x0 = [subjData.EF1.PxSubjMean(1);
                subjData.EF1.PySubjMean(1);
                subjData.EF1.VxSubjMean(1);
                subjData.EF1.VySubjMean(1);
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
EF1ChanParams.trialLength = length(subjDataChan.EF1.PxSubjMean);
EF1ChanParams.x0 = [subjDataChan.EF1.PxSubjMean(1);
                    subjDataChan.EF1.PySubjMean(1);
                    subjDataChan.EF1.VxSubjMean(1);
                    subjDataChan.EF1.VySubjMean(1);
                    subjDataChan.EF1.FxSubjMean(1);
                    subjDataChan.EF1.FySubjMean(1);
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
LF2Params.trialLength = length(subjData.LF2.PxSubjMean);
forceInit_LF2 = Params.Mass*[subjData.LF2.AxSubjMean(1); 
                            subjData.LF2.AySubjMean(1)] - ...
                LF2Params.curlGain*[0 1;-1 0]*[subjData.LF2.VxSubjMean(1);
                                               subjData.LF2.VySubjMean(1)];
LF2Params.x0 = [subjData.LF2.PxSubjMean(1);
            subjData.LF2.PySubjMean(1);
            subjData.LF2.VxSubjMean(1);
            subjData.LF2.VySubjMean(1);
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
LF2ChanParams.trialLength = length(subjDataChan.LF2.PxSubjMean);
LF2ChanParams.x0 = [subjDataChan.LF2.PxSubjMean(1);
                    subjDataChan.LF2.PySubjMean(1);
                    subjDataChan.LF2.VxSubjMean(1);
                    subjDataChan.LF2.VySubjMean(1);
                    subjDataChan.LF2.FxSubjMean(1);
                    subjDataChan.LF2.FySubjMean(1);
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
EN2Params.trialLength = length(subjData.EN2.PxSubjMean);
forceInit_EN2 = Params.Mass*[subjData.EN2.AxSubjMean(1); 
                            subjData.EN2.AySubjMean(1)] - ...
                EN2Params.curlGain*[0 1;-1 0]*[subjData.EN2.VxSubjMean(1);
                                               subjData.EN2.VySubjMean(1)];
EN2Params.x0 = [subjData.EN2.PxSubjMean(1);
            subjData.EN2.PySubjMean(1);
            subjData.EN2.VxSubjMean(1);
            subjData.EN2.VySubjMean(1);
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
EN2ChanParams.trialLength = length(subjDataChan.EN2.PxSubjMean);
EN2ChanParams.x0 = [subjDataChan.EN2.PxSubjMean(1);
                    subjDataChan.EN2.PySubjMean(1);
                    subjDataChan.EN2.VxSubjMean(1);
                    subjDataChan.EN2.VySubjMean(1);
                    subjDataChan.EN2.FxSubjMean(1);
                    subjDataChan.EN2.FySubjMean(1);
                    0;
                    0.1;0;0];
SOL_EN2Chan = generateReachTrajectory_Channel(EN2ChanParams);
ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_EN2Chan.x(:,1)';SOL_EN2Chan.x(:,2)';SOL_EN2Chan.x(:,3)';SOL_EN2Chan.x(:,4)'];
SOL_EN2Chan.Fx = ForceMat(1,:)';
SOL_EN2Chan.Fy = ForceMat(2,:)';
SOL_EN2Chan.AdaptInd= fitgain(SOL_EN2Chan.Fx,SOL_EN2Chan.x(:,4)); 
[~,t_ind] = max(abs(SOL_EN2Chan.Fx));
SOL_EN2Chan.FxMax = SOL_EN2Chan.Fx(t_ind);

%LN1
[~,mp_ind_model_LN1] = max(abs(SOL_LN1.x(:,1)));
[~,mp_ind_data_LN1] = max(abs(subjData.LN1.PxSubjMean));
[~,mp_ind_chandata_LN1] = max(abs(subjDataChan.LN1.FxSubjMean));

% RMS_LN1 = (((subjData.LN1.PxSubjMean(end) - SOL_LN1.x(end,1))./subjData.LN1.PxSubjSE(end))^2)/ci_scale + ...
%       (((subjData.LN1.PySubjMean(end) - SOL_LN1.x(end,2))./subjData.LN1.PySubjSE(end))^2)/ci_scale + ... 
%       10*(((subjData.LN1.PxSubjMean(mp_ind_data_LN1) - SOL_LN1.x(mp_ind_model_LN1,1))./subjData.LN1.PxSubjSE(mp_ind_data_LN1))^2)/ci_scale + ...
%       8*(((subjData.LN1.PySubjMean(mp_ind_data_LN1) - SOL_LN1.x(mp_ind_model_LN1,2))./subjData.LN1.PxSubjSE(mp_ind_data_LN1))^2)/ci_scale + ...
%       28*(((subjDataChan.LN1.FxSubjMean(mp_ind_chandata_LN1) - SOL_LN1Chan.FxMax)./subjDataChan.LN1.FxSubjSE(mp_ind_chandata_LN1))^2)/ci_scale + ...
%       23*(((subjDataChan.LN1.estGainMean - SOL_LN1Chan.AdaptInd)./subjDataChan.LN1.estGainSE)^2)/ci_scale;

RMS_LN1 = (((subjData.LN1.PxSubjMean(end) - SOL_LN1.x(end,1))./subjData.LN1.PxSubjSE(end))^2)/ci_scale^2 + ...
      (((subjData.LN1.PySubjMean(end) - SOL_LN1.x(end,2))./subjData.LN1.PySubjSE(end))^2)/ci_scale^2 + ... 
      10*(((subjData.LN1.PxSubjMean(mp_ind_data_LN1) - SOL_LN1.x(mp_ind_model_LN1,1))./subjData.LN1.PxSubjSE(mp_ind_data_LN1))^2)/ci_scale^2 + ...
      8*(((subjData.LN1.PySubjMean(mp_ind_data_LN1) - SOL_LN1.x(mp_ind_model_LN1,2))./subjData.LN1.PySubjSE(mp_ind_data_LN1))^2)/ci_scale^2 + ...
      28*(((subjDataChan.LN1.FxSubjMean(mp_ind_chandata_LN1) - SOL_LN1Chan.FxMax)./subjDataChan.LN1.FxSubjSE(mp_ind_chandata_LN1))^2)/ci_scale^2 + ...
      23*(((subjDataChan.LN1.estGainMean - SOL_LN1Chan.AdaptInd)./subjDataChan.LN1.estGainSE)^2)/ci_scale^2;


%EF1
[~,mp_ind_model_EF1] = max(abs(SOL_EF1.x(:,1)));
[~,mp_ind_data_EF1] = max(abs(subjData.EF1.PxSubjMean));
[~,mp_ind_chandata_EF1] = max(abs(subjDataChan.EF1.FxSubjMean));

RMS_EF1 = (((subjData.EF1.PxSubjMean(end) - SOL_EF1.x(end,1))./subjData.EF1.PxSubjSE(end))^2)/ci_scale^2 + ...
      (((subjData.EF1.PySubjMean(end) - SOL_EF1.x(end,2))./subjData.EF1.PySubjSE(end))^2)/ci_scale^2 + ... 
      10*(((subjData.EF1.PxSubjMean(mp_ind_data_EF1) - SOL_EF1.x(mp_ind_model_EF1,1))./subjData.EF1.PxSubjSE(mp_ind_data_EF1))^2)/ci_scale^2 + ...
      8*(((subjData.EF1.PySubjMean(mp_ind_data_EF1) - SOL_EF1.x(mp_ind_model_EF1,2))./subjData.EF1.PySubjSE(mp_ind_data_EF1))^2)/ci_scale^2 + ...
      28*(((subjDataChan.EF1.FxSubjMean(mp_ind_chandata_EF1) - SOL_EF1Chan.FxMax)./subjDataChan.EF1.FxSubjSE(mp_ind_chandata_EF1))^2)/ci_scale^2 + ...
      23*(((subjDataChan.EF1.estGainMean - SOL_EF1Chan.AdaptInd)./subjDataChan.EF1.estGainSE)^2)/ci_scale^2;

%LF2
[~,mp_ind_model_LF2] = max(abs(SOL_LF2.x(:,1)));
[~,mp_ind_data_LF2] = max(abs(subjData.LF2.PxSubjMean));
[~,mp_ind_chandata_LF2] = max(abs(subjDataChan.LF2.FxSubjMean));

% Squared, with AdaptInd
RMS_LF2 = (((subjData.LF2.PxSubjMean(end) - SOL_LF2.x(end,1))./subjData.LF2.PxSubjSE(end))^2)/ci_scale^2 + ...
      (((subjData.LF2.PySubjMean(end) - SOL_LF2.x(end,2))./subjData.LF2.PySubjSE(end))^2)/ci_scale^2 + ... 
      10*(((subjData.LF2.PxSubjMean(mp_ind_data_LF2) - SOL_LF2.x(mp_ind_model_LF2,1))./subjData.LF2.PxSubjSE(mp_ind_data_LF2))^2)/ci_scale^2 + ...
      8*(((subjData.LF2.PySubjMean(mp_ind_data_LF2) - SOL_LF2.x(mp_ind_model_LF2,2))./subjData.LF2.PySubjSE(mp_ind_data_LF2))^2)/ci_scale^2 + ...
      28*(((subjDataChan.LF2.FxSubjMean(mp_ind_chandata_LF2) - SOL_LF2Chan.FxMax)./subjDataChan.LF2.FxSubjSE(mp_ind_chandata_LF2))^2)/ci_scale^2 + ...
      23*(((subjDataChan.LF2.estGainMean - SOL_LF2Chan.AdaptInd)./subjDataChan.LF2.estGainSE)^2)/ci_scale^2;

%EN2
[~,mp_ind_model_EN2] = max(abs(SOL_EN2.x(:,1)));
[~,mp_ind_data_EN2] = max(abs(subjData.EN2.PxSubjMean));
[~,mp_ind_chandata_EN2] = max(abs(subjDataChan.EN2.FxSubjMean));

RMS_EN2 = (((subjData.EN2.PxSubjMean(end) - SOL_EN2.x(end,1))./subjData.EN2.PxSubjSE(end))^2)/ci_scale^2 + ...
      (((subjData.EN2.PySubjMean(end) - SOL_EN2.x(end,2))./subjData.EN2.PySubjSE(end))^2)/ci_scale^2 + ...
      10*(((subjData.EN2.PxSubjMean(mp_ind_data_EN2) - SOL_EN2.x(mp_ind_model_EN2,1))./subjData.EN2.PxSubjSE(mp_ind_data_EN2))^2)/ci_scale^2 + ...
      8*(((subjData.EN2.PySubjMean(mp_ind_data_EN2) - SOL_EN2.x(mp_ind_model_EN2,2))./subjData.EN2.PySubjSE(mp_ind_data_EN2))^2)/ci_scale^2 + ...
      28*(((subjDataChan.EN2.FxSubjMean(mp_ind_chandata_EN2) - SOL_EN2Chan.FxMax)./subjDataChan.EN2.FxSubjSE(mp_ind_chandata_EN2))^2)/ci_scale^2 + ...
      23*(((subjDataChan.EN2.estGainMean - SOL_EN2Chan.AdaptInd)./subjDataChan.EN2.estGainSE)^2)/ci_scale^2;

% Add factor so that if it is outside 95% confidence interval
% penalty = 0;
% if (((subjData.LF2.PxSubjMean(mp_ind_data_LF2) - SOL_LF2.x(mp_ind_model_LF2,1))./subjData.LF2.PxSubjSE(mp_ind_data_LF2))^2)/ci_scale > 1
%     penalty = penalty + 1;
% end
% if (((subjDataChan.LF2.FxSubjMean(mp_ind_chandata_LF2) - SOL_LF2Chan.FxMax)./subjDataChan.LF2.FxSubjSE(mp_ind_chandata_LF2))^2)/ci_scale > 1
%     penalty = penalty + 1;
% end
% if (((subjDataChan.LF2.estGainMean - SOL_LF2Chan.AdaptInd)./subjDataChan.LF2.estGainSE)^2)/ci_scale > 1
%     penalty = penalty + 1;
% end
% if (((subjData.EN2.PxSubjMean(mp_ind_data_EN2) - SOL_EN2.x(mp_ind_model_EN2,1))./subjData.EN2.PxSubjSE(mp_ind_data_EN2))^2)/ci_scale > 1
%     penalty = penalty + 1;
% end
% if (((subjDataChan.EN2.FxSubjMean(mp_ind_chandata_EN2) - SOL_EN2Chan.FxMax)./subjDataChan.EN2.FxSubjSE(mp_ind_chandata_EN2))^2)/ci_scale > 1
%     penalty = penalty + 1;
% end
% if (((subjDataChan.EN2.estGainMean - SOL_EN2Chan.AdaptInd)./subjDataChan.EN2.estGainSE)^2)/ci_scale > 1
%     penalty = penalty + 1;
% end
penalty = 0;
if ((subjData.LF2.PxSubjMean(mp_ind_data_LF2) - SOL_LF2.x(mp_ind_model_LF2,1))/subjData.LF2.PxSubjSE(mp_ind_data_LF2))^2 > ci_scale^2
    penalty = penalty + 1;
end
if ((subjDataChan.LF2.FxSubjMean(mp_ind_chandata_LF2) - SOL_LF2Chan.FxMax)/subjDataChan.LF2.FxSubjSE(mp_ind_chandata_LF2))^2 > ci_scale^2
    penalty = penalty + 1;
end
if ((subjDataChan.LF2.estGainMean - SOL_LF2Chan.AdaptInd)/subjDataChan.LF2.estGainSE)^2 > ci_scale^2
    penalty = penalty + 1;
end
if ((subjData.EN2.PxSubjMean(mp_ind_data_EN2) - SOL_EN2.x(mp_ind_model_EN2,1))/subjData.EN2.PxSubjSE(mp_ind_data_EN2))^2 > ci_scale^2
    penalty = penalty + 1;
end
if ((subjDataChan.EN2.FxSubjMean(mp_ind_chandata_EN2) - SOL_EN2Chan.FxMax)/subjDataChan.EN2.FxSubjSE(mp_ind_chandata_EN2))^2 > ci_scale^2
    penalty = penalty + 1;
end
if ((subjDataChan.EN2.estGainMean - SOL_EN2Chan.AdaptInd)/subjDataChan.EN2.estGainSE)^2 > ci_scale^2
    penalty = penalty + 1;
end

if penalty == 0
    % LF2 and EN2 meets z-score constraints
    X
    pause
end

% Combine Phases  
RMS = RMS_LN1 + RMS_EF1 + 5*RMS_LF2 + 5*RMS_EN2;