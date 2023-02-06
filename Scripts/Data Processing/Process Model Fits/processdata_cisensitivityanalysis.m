%% Process Data for Fig 4 and 6 -- CI sensitivity

ci = [95,90,80,70,60];
ages = {'Y','O'};
alpha = 0.4:0.05:1.2;


for ii = 1:length(ages)
    clearvars -except ages ii YoungAlphaSweep alpha ci CISweep Alphas
    for kk = 1:length(ci)
        if ages{ii} == 'Y'
            data_dir = ['../../../Data/Model Fits/Model Sensitivity Fits (Fig 5 and 6)/All Fits/Younger/',num2str(ci(kk)),'CI']; 
        elseif ages{ii} == 'O'
            data_dir = ['../../../Data/Model Fits/Model Sensitivity Fits (Fig 5 and 6)/All Fits/Older/',num2str(ci(kk)),'CI']; 
        else 
            error('Improper Input')
        end
        goodalphas = [];
        YoungAlphaSweep = [];
        for nn = 1:length(alpha)
            d = dir([data_dir,'/',ages{ii},'-1-alpha',num2str(alpha(nn)),'-*']);
            if ~isempty(d) && isempty(strfind(d(1).name,'badfit'))
                load([data_dir,'/',d(1).name],'Young')
                YoungAlphaSweep = [YoungAlphaSweep;
                                   Young];
                goodalphas = [goodalphas,alpha(nn)];
                load([data_dir,'/',d(1).name],'yData','yDataChan')

            else
                disp(['Sol ',num2str(alpha(nn)),'not found.'])
            end
        end
        CISweep(kk).YoungAlphaSweep = YoungAlphaSweep;
        CISweep(kk).goodalphas = goodalphas;
    end
    if ages{ii} == 'Y'
        save('../../../Data/Model Fits/Model Sensitivity Fits (Fig 5 and 6)/CISensitivity-Younger.mat','yData','yDataChan','CISweep')
    else
        save('../../../Data/Model Fits/Model Sensitivity Fits (Fig 5 and 6)/CISensitivity-Older.mat','yData','yDataChan','CISweep')
    end
end