function RMS = qualityOfModelFit_LF2EN2Only_individual(X,subjData,subjDataChan,FitIndices,Params,subjNum)
%% Calculate Quality of Model Fit
% This function is what is minimized 
% (smaller number ~ better fitting model)
%
% This particular version only considers the quality of fit for Late
% Learning (LF2) and Early Washout (EN2). 
%
% Data from both Channel Trials and Regular trials are used
% 
% Quality of Fit looks at how closely the position of maximum perpendicular
% error matches, the end point, the maximum perpendicular force, and
% adaptation index. 
%
% EN2 and LF2 are equally weighted.
% 
% This looks an individual subject fits without confidence bounds


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
elseif FitIndices(11)
    Params.estGainEN2 = Params.estGainLF2;
else
    Params.estGainEN2 = Params.estGain;
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

% LF2
LF2Params = Params;
if FitIndices(11)
    LF2Params.estGain = Params.estGainLF2;
else
    LF2Params.estGain = Params.estGain;    
end
LF2Params.curlGain = -20;


% LF2Params.trialLength = find(~isnan(subjData.LF2.PxSubjNorm(:,subjNum)),1,'last');
% LF2Params.trialLength = subjData.LF2.SubjTrLength(subjNum);
LF2Params.trialLength = floor(subjData.LF2.SubjMeanTrLength(subjNum));

% subj_ind_LF2 = find(subjData.LF2.SubjTrLength(subjNum) == nonzeros(subjData.LF2.SubjTrLength),1,'first');
% subj_ind_LF2 = find(subjData.LF2.SubjMeanTrLength(subjNum) == nonzeros(subjData.LF2.SubjMeanTrLength),1,'first');
% subj_ind_LF2 = find(subjData.LF2.SubjMeanTrLength(subjNum) == subjData.LF2.SubjMeanTrLength(~isnan(subjData.LF2.SubjMeanTrLength)),1,'first');
if isnan(subjData.LF2.SubjMeanTrLength(subjNum))
    subj_ind_LF2 = NaN;
else
    subj_ind_LF2 = sum(~isnan(subjData.LF2.SubjMeanTrLength(1:subjNum)));
end

forceInit_LF2 = Params.Mass*[subjData.LF2.AxSubjNorm(1,subj_ind_LF2); 
                            subjData.LF2.AySubjNorm(1,subj_ind_LF2)] - ...
                LF2Params.curlGain*[0 1;-1 0]*[subjData.LF2.VxSubjNorm(1,subj_ind_LF2);
                                               subjData.LF2.VySubjNorm(1,subj_ind_LF2)];
LF2Params.x0 = [subjData.LF2.PxSubjNorm(1,subj_ind_LF2);
            subjData.LF2.PySubjNorm(1,subj_ind_LF2);
            subjData.LF2.VxSubjNorm(1,subj_ind_LF2);
            subjData.LF2.VySubjNorm(1,subj_ind_LF2);
            forceInit_LF2(1);forceInit_LF2(2);
            0;
            0.1;
            0;0];
SOL_LF2 = generateReachTrajectory_CurlOrBaseline(LF2Params);
AccMat = inv(Params.Mass)*([SOL_LF2.x(:,5)';SOL_LF2.x(:,6)'] + LF2Params.curlGain*[0 1;-1 0]*[SOL_LF2.x(:,3)';SOL_LF2.x(:,4)']);
SOL_LF2.Ax = AccMat(1,:)';
SOL_LF2.Ay = AccMat(2,:)';


if isnan(subjDataChan.LF2.SubjMeanTrLength(subjNum))
    subj_ind_LF2_Chan = NaN;
else
    subj_ind_LF2_Chan = sum(~isnan(subjDataChan.LF2.SubjMeanTrLength(1:subjNum)));
end

% if subjDataChan.LF2.SubjTrLength(subjNum) ~= 0
if ~isnan(subjDataChan.LF2.SubjMeanTrLength(subjNum))
    LF2_Chan_OK = 1;
%     subj_ind_LF2_Chan = find(subjDataChan.LF2.SubjTrLength(subjNum) == nonzeros(subjDataChan.LF2.SubjTrLength),1,'first');
%     subj_ind_LF2_Chan = find(subjDataChan.LF2.SubjMeanTrLength(subjNum) == nonzeros(subjDataChan.LF2.SubjMeanTrLength),1,'first');
%     subj_ind_LF2_Chan = find(subjDataChan.LF2.SubjMeanTrLength(subjNum) == subjDataChan.LF2.SubjMeanTrLength(~isnan(subjDataChan.LF2.SubjMeanTrLength)),1,'first');
    % LF2 CHANNEL
    LF2ChanParams = Params;
    if FitIndices(11)
        LF2ChanParams.estGain = Params.estGainLF2;
    else
        LF2ChanParams.estGain = Params.estGain;
    end
    LF2ChanParams.curlGain = -20;
%     LF2ChanParams.trialLength = find(~isnan(subjDataChan.LF2.PxSubjNorm(:,subjNum)),1,'last');
    LF2ChanParams.trialLength = find(~isnan(subjDataChan.LF2.PxSubjNorm(:,subj_ind_LF2_Chan)),1,'last');
    LF2ChanParams.x0 = [subjDataChan.LF2.PxSubjNorm(1,subj_ind_LF2_Chan);
                        subjDataChan.LF2.PySubjNorm(1,subj_ind_LF2_Chan);
                        subjDataChan.LF2.VxSubjNorm(1,subj_ind_LF2_Chan);
                        subjDataChan.LF2.VySubjNorm(1,subj_ind_LF2_Chan);
                        subjDataChan.LF2.FxSubjNorm(1,subj_ind_LF2_Chan);
                        subjDataChan.LF2.FySubjNorm(1,subj_ind_LF2_Chan);
                        0;
                        0.1;0;0];
    SOL_LF2Chan = generateReachTrajectory_Channel(LF2ChanParams);
    ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_LF2Chan.x(:,1)';SOL_LF2Chan.x(:,2)';SOL_LF2Chan.x(:,3)';SOL_LF2Chan.x(:,4)'];
    SOL_LF2Chan.Fx = ForceMat(1,:)';
    SOL_LF2Chan.Fy = ForceMat(2,:)';
    SOL_LF2Chan.AdaptInd= fitgain(SOL_LF2Chan.Fx,SOL_LF2Chan.x(:,4)); 
    [~,t_ind] = max(abs(SOL_LF2Chan.Fx));
    SOL_LF2Chan.FxMax = SOL_LF2Chan.Fx(t_ind);
else
    LF2_Chan_OK = 0;
end

% EN2
EN2Params = Params;
if FitIndices(12)
    EN2Params.estGain = Params.estGainEN2;
elseif FitIndices(11)
    EN2Params.estGain = Params.estGainLF2;
else
    EN2Params.estGain = subjData.EN2.estGainMean;
end
EN2Params.curlGain = 0;
% EN2Params.trialLength = find(~isnan(subjData.EN2.PxSubjNorm(:,subjNum)),1,'last');
% EN2Params.trialLength = subjData.EN2.SubjTrLength(subjNum);
EN2Params.trialLength = floor(subjData.EN2.SubjMeanTrLength(subjNum));

% subj_ind_EN2 = find(subjData.EN2.SubjTrLength(subjNum) == nonzeros(subjData.EN2.SubjTrLength),1,'first');
% subj_ind_EN2 = find(subjData.EN2.SubjMeanTrLength(subjNum) == nonzeros(subjData.EN2.SubjMeanTrLength),1,'first');
% subj_ind_EN2 = find(subjData.EN2.SubjMeanTrLength(subjNum) == subjData.EN2.SubjMeanTrLength(~isnan(subjData.EN2.SubjMeanTrLength)),1,'first');
if isnan(subjData.EN2.SubjMeanTrLength(subjNum))
    subj_ind_EN2 = NaN;
else
    subj_ind_EN2 = sum(~isnan(subjData.EN2.SubjMeanTrLength(1:subjNum)));
end

forceInit_EN2 = Params.Mass*[subjData.EN2.AxSubjNorm(1,subj_ind_EN2); 
                            subjData.EN2.AySubjNorm(1,subj_ind_EN2)] - ...
                EN2Params.curlGain*[0 1;-1 0]*[subjData.EN2.VxSubjNorm(1,subj_ind_EN2);
                                               subjData.EN2.VySubjNorm(1,subj_ind_EN2)];
EN2Params.x0 = [subjData.EN2.PxSubjNorm(1,subj_ind_EN2);
            subjData.EN2.PySubjNorm(1,subj_ind_EN2);
            subjData.EN2.VxSubjNorm(1,subj_ind_EN2);
            subjData.EN2.VySubjNorm(1,subj_ind_EN2);
            forceInit_EN2(1);forceInit_EN2(2);
            0;
            0.1;
            0;0];
SOL_EN2 = generateReachTrajectory_CurlOrBaseline(EN2Params);
AccMat = inv(Params.Mass)*([SOL_EN2.x(:,5)';SOL_EN2.x(:,6)'] + EN2Params.curlGain*[0 1;-1 0]*[SOL_EN2.x(:,3)';SOL_EN2.x(:,4)']);
SOL_EN2.Ax = AccMat(1,:)';
SOL_EN2.Ay = AccMat(2,:)';


if isnan(subjDataChan.EN2.SubjMeanTrLength(subjNum))
    subj_ind_EN2_Chan = NaN;
else
    subj_ind_EN2_Chan = sum(~isnan(subjDataChan.EN2.SubjMeanTrLength(1:subjNum)));
end
% if subjDataChan.EN2.SubjTrLength(subjNum) ~= 0
if ~isnan(subjDataChan.EN2.SubjMeanTrLength(subjNum))
%     subj_ind_EN2_Chan = find(subjDataChan.EN2.SubjTrLength(subjNum) == nonzeros(subjDataChan.EN2.SubjTrLength),1,'first');
%     subj_ind_EN2_Chan = find(subjDataChan.EN2.SubjMeanTrLength(subjNum) == nonzeros(subjDataChan.EN2.SubjMeanTrLength),1,'first');
%     subj_ind_EN2_Chan = find(subjDataChan.EN2.SubjMeanTrLength(subjNum) == subjDataChan.EN2.SubjMeanTrLength(~isnan(subjDataChan.EN2.SubjMeanTrLength)),1,'first');



    % EN2 CHANNEL
    EN2ChanParams = Params;
    if FitIndices(12)
        EN2ChanParams.estGain = Params.estGainEN2;
    elseif FitIndices(11)
        EN2ChanParams.estGain = Params.estGainLF2;
    else
        EN2ChanParams.estGain = subjDataChan.EN2.estGainMean;
    end
    EN2ChanParams.curlGain = 0;
    EN2ChanParams.trialLength = find(~isnan(subjDataChan.EN2.PxSubjNorm(:,subj_ind_EN2_Chan)),1,'last');
    EN2ChanParams.x0 = [subjDataChan.EN2.PxSubjNorm(1,subj_ind_EN2_Chan);
                        subjDataChan.EN2.PySubjNorm(1,subj_ind_EN2_Chan);
                        subjDataChan.EN2.VxSubjNorm(1,subj_ind_EN2_Chan);
                        subjDataChan.EN2.VySubjNorm(1,subj_ind_EN2_Chan);
                        subjDataChan.EN2.FxSubjNorm(1,subj_ind_EN2_Chan);
                        subjDataChan.EN2.FySubjNorm(1,subj_ind_EN2_Chan);
                        0;
                        0.1;0;0];
    SOL_EN2Chan = generateReachTrajectory_Channel(EN2ChanParams);
    ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_EN2Chan.x(:,1)';SOL_EN2Chan.x(:,2)';SOL_EN2Chan.x(:,3)';SOL_EN2Chan.x(:,4)'];
    SOL_EN2Chan.Fx = ForceMat(1,:)';
    SOL_EN2Chan.Fy = ForceMat(2,:)';
    SOL_EN2Chan.AdaptInd= fitgain(SOL_EN2Chan.Fx,SOL_EN2Chan.x(:,4)); 
    [~,t_ind] = max(abs(SOL_EN2Chan.Fx));
    SOL_EN2Chan.FxMax = SOL_EN2Chan.Fx(t_ind);
else
    EN2_Chan_OK = 0;
end

%% Calculate Objective Function

%LF2
[~,mpx_ind] = max(abs(subjData.LF2.PxSubjNorm(:,subj_ind_LF2)));
[~,mpx_ind_mdl] = max(abs(SOL_LF2.x(:,1)));
% [~,mp_ind_chandata] = max(abs(subjDataChan.LF2.FxSubjNorm(:,subj_ind_LF2_Chan)));
endpt = find(~isnan(subjData.LF2.PxSubjNorm(:,subj_ind_LF2)),1,'last');


% if isnan(subjDataChan.LF2.FxMax(subjNum)) || LF2_Chan_OK == 0
    RMS_LF2 = (subjData.LF2.PxSubjNorm(endpt,subj_ind_LF2) - SOL_LF2.x(end,1))^2 + ... % End Px
          (subjData.LF2.PySubjNorm(endpt,subj_ind_LF2) - SOL_LF2.x(end,2))^2 + ... % End Py
          (subjData.LF2.PxSubjNorm(mpx_ind,subj_ind_LF2) - SOL_LF2.x(mpx_ind_mdl,1))^2 + ... % Max Px
          (subjData.LF2.PySubjNorm(mpx_ind,subj_ind_LF2) - SOL_LF2.x(mpx_ind_mdl,2))^2; % Py of Max Px
% else 
%     RMS_LF2 = (subjData.LF2.PxSubjNorm(endpt,subj_ind_LF2) - SOL_LF2.x(end,1))^2 + ... % End Px
%           (subjData.LF2.PySubjNorm(endpt,subj_ind_LF2) - SOL_LF2.x(end,2))^2 + ... % End Py
%           15*(subjData.LF2.PxSubjNorm(mpx_ind,subj_ind_LF2) - SOL_LF2.x(mpx_ind_mdl,1))^2 + ... % Max Px
%           8*(subjData.LF2.PySubjNorm(mpx_ind,subj_ind_LF2) - SOL_LF2.x(mpx_ind_mdl,2))^2 + ... % Py of Max Px
%           1e-4*(subjDataChan.LF2.FxSubjNorm(mp_ind_chandata,subj_ind_LF2_Chan) - SOL_LF2Chan.FxMax)^2 + ... % Max Fx
%           1e-4*(subjDataChan.LF2.estGain(subjNum) - SOL_LF2Chan.AdaptInd)^2; % Adapt Ind
% end

%EN2
[~,mpx_ind] = max(abs(subjData.EN2.PxSubjNorm(:,subj_ind_EN2)));
[~,mpx_ind_mdl] = max(abs(SOL_EN2.x(:,1)));
% [~,mp_ind_chandata] = max(abs(subjDataChan.EN2.FxSubjNorm(:,subj_ind_EN2_Chan)));
endpt = find(~isnan(subjData.EN2.PxSubjNorm(:,subj_ind_EN2)),1,'last');

% if isnan(subjDataChan.EN2.FxMax(subjNum)) || EN2_Chan_OK == 0
    RMS_EN2 = (subjData.EN2.PxSubjNorm(endpt,subj_ind_EN2) - SOL_EN2.x(end,1))^2 + ... % End Px
          (subjData.EN2.PySubjNorm(endpt,subj_ind_EN2) - SOL_EN2.x(end,2))^2 + ... % End Py
          (subjData.EN2.PxSubjNorm(mpx_ind,subj_ind_EN2) - SOL_EN2.x(mpx_ind_mdl,1))^2 + ... % Max Px
          (subjData.EN2.PySubjNorm(mpx_ind,subj_ind_EN2) - SOL_EN2.x(mpx_ind_mdl,2))^2; % Py of Max Px
% else
%     RMS_EN2 = (subjData.EN2.PxSubjNorm(endpt,subj_ind_EN2) - SOL_EN2.x(end,1))^2 + ... % End Px
%           (subjData.EN2.PySubjNorm(endpt,subj_ind_EN2) - SOL_EN2.x(end,2))^2 + ... % End Py
%           15*(subjData.EN2.PxSubjNorm(mpx_ind,subj_ind_EN2) - SOL_EN2.x(mpx_ind_mdl,1))^2 + ... % Max Px
%           8*(subjData.EN2.PySubjNorm(mpx_ind,subj_ind_EN2) - SOL_EN2.x(mpx_ind_mdl,2))^2 + ... % Py of Max Px
%           1e-4*(subjDataChan.EN2.FxSubjNorm(mp_ind_chandata,subj_ind_EN2_Chan) - SOL_EN2Chan.FxMax)^2 + ... % Max Fx
%           1e-4*(subjDataChan.EN2.estGain(subjNum) - SOL_EN2Chan.AdaptInd)^2; % Adapt Ind
% end


RMS = RMS_LF2 + RMS_EN2;
if isnan(RMS)
    pause
end
