%% Produce Processed Data for Both Older and Younger Adults

% Younger
yData = processRawCurlTrialData('Y',2);
yDataChan = processRawChannelTrialData('Y',2);
save('../../../Data/Trajectory Data/Processed Data/Younger-AveragedTrajectoryData.mat','yData','yDataChan')
% save('../../Data/Trajectory Data/Processed Data/Younger-AveragedTrajectoryData-Individual.mat','yData','yDataChan')

% Older
oData = processRawCurlTrialData('O',2);
oDataChan = processRawChannelTrialData('O',2);
clear yData yDataChan yBootData
yData = oData;
yDataChan = oDataChan;
save('../../../Data/Trajectory Data/Processed Data/Older-AveragedTrajectoryData-10tr.mat','yData','yDataChan')
% save('../../../Data/Trajectory Data/Processed Data/Older-AveragedTrajectoryData.mat','yData','yDataChan')
% save('../../Data/Trajectory Data/Processed Data/Older-AveragedTrajectoryData-Individual.mat','yData','yDataChan','yBootData')