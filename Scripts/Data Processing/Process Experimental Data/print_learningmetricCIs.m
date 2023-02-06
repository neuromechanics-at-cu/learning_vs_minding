% Print all the error metrics

phaseNames = {'LN1','EF1','LF2','EN2'};
ages = {'Y','O'};

for jj = 1:length(ages)
    %% Load Data
    if ages{jj} == 'Y'
        load('../../../Data/Trajectory Data/Processed Data/Younger-AveragedTrajectoryData.mat')
    else
        load('../../../Data/Trajectory Data/Processed Data/Older-AveragedTrajectoryData.mat')
    end
    disp(' --- ')
    disp(ages{jj})

    for ii = 1:length(phaseNames)
        eval(['[~,mp_ind] = max(abs(yData.',phaseNames{ii},'.PxSubjMean));'])
        eval(['PxMean = yData.',phaseNames{ii},'.PxSubjMean(mp_ind);'])
        eval(['PxSE = yData.',phaseNames{ii},'.PxSubjSE(mp_ind);'])
        disp([phaseNames{ii},', Px: ',num2str(PxMean*100),' +/- ',num2str(PxSE*100)])
        eval(['[~,mp_ind] = max(abs(yDataChan.',phaseNames{ii},'.FxSubjMean));'])
        eval(['FxMean = yDataChan.',phaseNames{ii},'.FxSubjMean(mp_ind);'])
        eval(['FxSE = yDataChan.',phaseNames{ii},'.FxSubjSE(mp_ind);'])
        disp([phaseNames{ii},', Fx: ',num2str(FxMean),' +/- ',num2str(FxSE)])
        eval(['adaptIndMean = yDataChan.',phaseNames{ii},'.estGainMean;'])
        eval(['adaptIndSE = yDataChan.',phaseNames{ii},'.estGainSE;'])
        disp([phaseNames{ii},', EstGain: ',num2str(adaptIndMean),' +/- ',num2str(adaptIndSE)])
        disp('  ')
    end
    disp(' ----- ')
end

phaseNames = {'LN1','EF1','LF2','EN2'};
ages = {'Y','O'};

disp('Indiviual Processed Data')
for jj = 1:length(ages)
    %% Load Data
    if ages{jj} == 'Y'
        load('../../../Data/Trajectory Data/Processed Data/Younger-AveragedTrajectoryData-Individual.mat')
    else
        load('../../../Data/Trajectory Data/Processed Data/Older-AveragedTrajectoryData-Individual.mat')
    end
    disp(' --- ')
    disp(ages{jj})

    for ii = 1:length(phaseNames)
        eval(['[~,mp_ind] = max(abs(yData.',phaseNames{ii},'.PxSubjMean));'])
        eval(['PxMean = yData.',phaseNames{ii},'.PxSubjMean(mp_ind);'])
        eval(['PxSE = yData.',phaseNames{ii},'.PxSubjSE(mp_ind);'])
        disp([phaseNames{ii},', Px: ',num2str(PxMean*100),' +/- ',num2str(PxSE*100)])
        eval(['[~,mp_ind] = max(abs(yDataChan.',phaseNames{ii},'.FxSubjMean));'])
        eval(['FxMean = yDataChan.',phaseNames{ii},'.FxSubjMean(mp_ind);'])
        eval(['FxSE = yDataChan.',phaseNames{ii},'.FxSubjSE(mp_ind);'])
        disp([phaseNames{ii},', Fx: ',num2str(FxMean),' +/- ',num2str(FxSE)])
        eval(['adaptIndMean = yDataChan.',phaseNames{ii},'.estGainMean;'])
        eval(['adaptIndSE = yDataChan.',phaseNames{ii},'.estGainSE;'])
        disp([phaseNames{ii},', EstGain: ',num2str(adaptIndMean),' +/- ',num2str(adaptIndSE)])
        disp('  ')
    end
    disp(' ----- ')
end