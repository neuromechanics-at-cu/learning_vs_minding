%% Calculate Quality of Model Fit Using Trajectory Position
% This function is what is minimized 
% (smaller number ~ better fitting model)
%
% Data from Regular trials are used
% 
% Quality of Fit looks at how closely the position vs time over the entire
% trajectory matches.
%
% LN1, EF1, EN2, and LF2 are equally weighted.

function [RMS] = qualityOfModelFit_AllPhases_PosVel(X,subjData,subjDataChan,FitIndices,Params)

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

%LN1
LN1Params = Params;
if FitIndices(9)
    LN1Params.estGain = X(count);
    count = count + 1;
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
            subjData.LN1.PxSubjMean(end);
            subjData.LN1.PySubjMean(end);0;0];
SOL_LN1 = generateReachTrajectory_CurlOrBaseline(LN1Params);
AccMat = inv(Params.Mass)*([SOL_LN1.x(:,5)';SOL_LN1.x(:,6)'] + LN1Params.curlGain*[0 1;-1 0]*[SOL_LN1.x(:,3)';SOL_LN1.x(:,4)']);
SOL_LN1.Ax = AccMat(1,:)';
SOL_LN1.Ay = AccMat(2,:)';


% EF1
EF1Params = Params;
if FitIndices(10)
    EF1Params.estGain = X(count);
    count = count + 1;
else
    EF1Params.estGain = subjData.EF1.estGainMean;
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
            subjData.EF1.PxSubjMean(end);
            subjData.EF1.PySubjMean(end);0;0];
SOL_EF1 = generateReachTrajectory_CurlOrBaseline(EF1Params);
AccMat = inv(Params.Mass)*([SOL_EF1.x(:,5)';SOL_EF1.x(:,6)'] + EF1Params.curlGain*[0 1;-1 0]*[SOL_EF1.x(:,3)';SOL_EF1.x(:,4)']);
SOL_EF1.Ax = AccMat(1,:)';
SOL_EF1.Ay = AccMat(2,:)';


% LF2
LF2Params = Params;
if FitIndices(11)
    LF2Params.estGain = X(count);
    count = count + 1;
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
            subjData.LF2.PxSubjMean(end);
            subjData.LF2.PySubjMean(end);0;0];
SOL_LF2 = generateReachTrajectory_CurlOrBaseline(LF2Params);
AccMat = inv(Params.Mass)*([SOL_LF2.x(:,5)';SOL_LF2.x(:,6)'] + LF2Params.curlGain*[0 1;-1 0]*[SOL_LF2.x(:,3)';SOL_LF2.x(:,4)']);
SOL_LF2.Ax = AccMat(1,:)';
SOL_LF2.Ay = AccMat(2,:)';


% EN2
EN2Params = Params;
if FitIndices(12)
    EN2Params.estGain = X(count);
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
            subjData.EN2.PxSubjMean(end);
            subjData.EN2.PySubjMean(end);0;0];
SOL_EN2 = generateReachTrajectory_CurlOrBaseline(EN2Params);
AccMat = inv(Params.Mass)*([SOL_EN2.x(:,5)';SOL_EN2.x(:,6)'] + EN2Params.curlGain*[0 1;-1 0]*[SOL_EN2.x(:,3)';SOL_EN2.x(:,4)']);
SOL_EN2.Ax = AccMat(1,:)';
SOL_EN2.Ay = AccMat(2,:)';



% Calculate Error
% LN1
RMS_LN1 = sum(abs((subjData.LN1.PxSubjMean' - SOL_LN1.x(:,1))./subjData.LN1.PxSubjSE')) + ...
          sum(abs((subjData.LN1.PySubjMean' - SOL_LN1.x(:,2))./subjData.LN1.PySubjSE')) + ...
          sum(abs((subjData.LN1.VxSubjMean' - SOL_LN1.x(:,3))./subjData.LN1.VxSubjSE')) + ...
          sum(abs((subjData.LN1.VySubjMean' - SOL_LN1.x(:,4))./subjData.LN1.VySubjSE'));
RMS_LN1 = RMS_LN1/4/LN1Params.trialLength;

% EF1
RMS_EF1 = sum(abs((subjData.EF1.PxSubjMean' - SOL_EF1.x(:,1))./subjData.EF1.PxSubjSE')) + ...
          sum(abs((subjData.EF1.PySubjMean' - SOL_EF1.x(:,2))./subjData.EF1.PySubjSE')) + ...
          sum(abs((subjData.EF1.VxSubjMean' - SOL_EF1.x(:,3))./subjData.EF1.VxSubjSE')) + ...
          sum(abs((subjData.EF1.VySubjMean' - SOL_EF1.x(:,4))./subjData.EF1.VySubjSE'));
RMS_EF1 = RMS_EF1/4/EF1Params.trialLength;

% LF2 
RMS_LF2 = sum(abs((subjData.LF2.PxSubjMean' - SOL_LF2.x(:,1))./subjData.LF2.PxSubjSE')) + ...
          sum(abs((subjData.LF2.PySubjMean' - SOL_LF2.x(:,2))./subjData.LF2.PySubjSE')) + ...
          sum(abs((subjData.LF2.VxSubjMean' - SOL_LF2.x(:,3))./subjData.LF2.VxSubjSE')) + ...
          sum(abs((subjData.LF2.VySubjMean' - SOL_LF2.x(:,4))./subjData.LF2.VySubjSE'));
RMS_LF2 = RMS_LF2/4/LF2Params.trialLength;

% EN2
RMS_EN2 = sum(abs((subjData.EN2.PxSubjMean' - SOL_EN2.x(:,1))./subjData.EN2.PxSubjSE')) + ...
          sum(abs((subjData.EN2.PySubjMean' - SOL_EN2.x(:,2))./subjData.EN2.PySubjSE')) + ...
          sum(abs((subjData.EN2.VxSubjMean' - SOL_EN2.x(:,3))./subjData.EN2.VxSubjSE')) + ...
          sum(abs((subjData.EN2.VySubjMean' - SOL_EN2.x(:,4))./subjData.EN2.VySubjSE'));
RMS_EN2 = RMS_EN2/4/EN2Params.trialLength;

RMS = RMS_LN1 + RMS_EF1 + RMS_LF2 + RMS_EN2;
