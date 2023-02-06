%% Model Fits Across a Range of Learning

% Make a matrix
ages = ['Y','O'];
ci = [95,90,80,70,60];
alpha_range = 0.6:0.05:0.9;
Results_Matrix = zeros(length(ci),length(alpha_range));
Results_Matrix_Y = zeros(length(ci),length(alpha_range));
Results_Matrix_O = zeros(length(ci),length(alpha_range));

for nn = 1:length(ages)
    clearvars -except ages nn ci alpha_range Results_Matrix Results_Matrix_Y Results_Matrix_O

    %% Load Data
    if ages(nn) == 'Y'
        load('../../Data/Model Fits/Model Sensitivity Fits (Fig 5 and 6)/CISensitivity-Younger.mat')
        val = 2;

    elseif ages(nn) == 'O'
        load('../../Data/Model Fits/Model Sensitivity Fits (Fig 5 and 6)/CISensitivity-Older.mat')
        val = -1;
    else 
        error('Improper Input')
    end

    for mm = 1:length(ci)
        for kk = 1:length(CISweep(mm).goodalphas)
            alpha_ind = find(abs(alpha_range - CISweep(mm).goodalphas(kk))<4*eps); 
            Results_Matrix(mm,alpha_ind) = Results_Matrix(mm,alpha_ind) + val;
            if nn == 1
                Results_Matrix_Y(mm,alpha_ind) = Results_Matrix_Y(mm,alpha_ind) + 1;
            elseif nn == 2
                Results_Matrix_O(mm,alpha_ind) = Results_Matrix_O(mm,alpha_ind) + 1;
            end
        end
    end
end
Results_Matrix(Results_Matrix==0) = NaN;
Results_Matrix(Results_Matrix==-1) = 0;

figure(601)
hm = heatmap(Results_Matrix);
hm.XDisplayLabels = string(alpha_range);
hm.YDisplayLabels = string(ci);
hm.ColorbarVisible = 'off';
hm.MissingDataColor = [0.8,0.8,0.8];
hm.CellLabelColor = 'none';
hm.YLabel = 'Confidence Interval';
hm.XLabel = 'Proportion Learned';
hm.Title = 'Model Analysis: Confidence Interval Limit';
colorMap = [0,104,56;0,128,128;46,49,146]./255;
colormap(colorMap)

% Plot for each age group (for illustrator processing)
if 0
    figure(602)
    hm = heatmap(Results_Matrix_Y);
    hm.XDisplayLabels = string(alpha_range);
    hm.YDisplayLabels = string(ci);
    hm.ColorbarVisible = 'off';
    hm.GridVisible = 'off';
    hm.MissingDataColor = [0.8,0.8,0.8];
    hm.CellLabelColor = 'none';
    hm.YLabel = 'Confidence Interval';
    hm.XLabel = 'Proportion Learned';
    hm.Title = 'Model Analysis: Confidence Interval Limit - Young';
    colorMap = [255,255,255;0,104,56]./255;
    colormap(colorMap)
    
    figure(603)
    hm = heatmap(Results_Matrix_O);
    hm.XDisplayLabels = string(alpha_range);
    hm.YDisplayLabels = string(ci);
    hm.ColorbarVisible = 'off';
    hm.GridVisible = 'off';
    hm.MissingDataColor = [0.8,0.8,0.8];
    hm.CellLabelColor = 'none';
    hm.YLabel = 'Confidence Interval';
    hm.XLabel = 'Proportion Learned';
    hm.Title = 'Model Analysis: Confidence Interval Limit - Young';
    colorMap = [255,255,255;46,49,146]./255;
    colormap(colorMap)
end

