%% Form Data Files for Figure 4
% alpha: range of alphas
% yData, yDataChan: processed data
% Young: is this used? "best fit model" from Fig 3? Not used in the plots.
% YoungAlphaSweep

ages = {'Y','O'};
alpha = 0.4:0.05:1.2;


for ii = 1:length(ages)
    clearvars -except ages ii YoungAlphaSweep alpha

    if ages{ii} == 'Y'
        data_dir = '../../../Data/Model Fits/Model Sensitivity Fits (Fig 5 and 6)/All Fits/Younger/95CI'; 
    elseif ages{ii} == 'O'
        data_dir = '../../../Data/Model Fits/Model Sensitivity Fits (Fig 5 and 6)/All Fits/Older/95CI'; 
    else 
        error('Improper Input')
    end

    for nn = 1:length(alpha)
        d = dir([data_dir,'/',ages{ii},'-1-alpha',num2str(alpha(nn)),'*']);
        if ~isempty(d) 
            load([data_dir,'/',d(1).name],'Young')
            YoungAlphaSweep(nn) = Young;
        else
            disp(['Sol ',num2str(alpha(nn)),'not found.'])
        end
    end

    load([data_dir,'/',d(1).name],'yData','yDataChan')
    if ages{ii} == 'Y'
        save('../../../Data/Model Fits/Model Sensitivity Fits (Fig 5 and 6)/ModelSensitivity-Younger.mat','yData','yDataChan','YoungAlphaSweep','alpha','Young')
    else
        save('../../../Data/Model Fits/Model Sensitivity Fits (Fig 5 and 6)/ModelSensitivity-Older.mat','yData','yDataChan','YoungAlphaSweep','alpha','Young')
    end
end

