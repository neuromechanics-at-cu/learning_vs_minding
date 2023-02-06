%% Print Table of Search Function Validation Results

% Load Data
load('../../Data/Model Fits/Reach Model and Search Validation (Fig 3)/SearchFunctionValidationResults.mat')

% Print
% Header
fprintf('True Learning \t True Cost Ratio \t Vary Learning? \t Vary Costs? \t Num Traj \t Found Learning \t Found Costs\n');
% Define +/- sign
pmchar = char(177);

%% Display All Results
% for ii = 1:length(CR_strings)
%     for jj = 1:length(learn_strings)
%         for kk = 1:length(searchtype_strings)
%             for ll = 1:length(numtraj_strings)
%                 temp = All_Searches{ii,jj,kk,ll};
%                 fprintf(['%1.3f \t\t %1.3f \t\t\t %1d \t\t\t %1d \t\t %1d \t\t %1.3f',pmchar,'%1.3f \t\t %1.3f',pmchar,'%1.3f \n'], ...
%                     temp.true_estgain/-20,temp.true_costratio,...
%                         temp.vary_learning,temp.vary_costs,temp.one_or_two_traj,...
%                         temp.estgain_mean/-20,temp.estgain_std/20,...
%                         temp.costratio_mean,temp.costratio_std);
%             end
%         end
%     end
% end

%% Display Just Corner Case
for ii = 1:length(CR_strings)
    for jj = 1:length(learn_strings)
        for kk = 1:length(searchtype_strings)
            for ll = 1:length(numtraj_strings)
                temp = All_Searches{ii,jj,kk,ll};
                if temp.vary_costs && temp.vary_learning && (temp.one_or_two_traj == 2)
                    fprintf(['%1.3f \t\t %1.3f \t\t\t %1d \t\t\t %1d \t\t %1d \t\t %1.3f',pmchar,'%1.3f \t\t %1.3f',pmchar,'%1.3f \n'], ...
                        temp.true_estgain/-20,temp.true_costratio,...
                        temp.vary_learning,temp.vary_costs,temp.one_or_two_traj,...
                        temp.estgain_mean/-20,temp.estgain_std/20,...
                        temp.costratio_mean,temp.costratio_std);
                    % Save These For Later
                    if temp.true_costratio < -3 && temp.true_estgain/-20 < 0.8
                        Temp_1 = temp;
                    elseif temp.true_costratio < -3 && temp.true_estgain/-20 > 0.8
                        Temp_2 = temp;
                    end
                end
            end
        end
    end
end

%% Plot Example Fits
% Choosing 60% and 100% for Low Cost Ratio

for nn = 1:2
    if nn == 1
        temp = Temp_1;
    else
        temp = Temp_2;
    end
    figure(315+nn)
    for kk = 1:10
        subplot(1,2,1)
        hold on
        plot(temp.AlphaSweep{kk}.SOL_LF2.x(:,1),temp.AlphaSweep{kk}.SOL_LF2.x(:,2),'color',[0.65,0.65,0.65])
        subplot(1,2,2)
        hold on
        plot(temp.AlphaSweep{kk}.SOL_EN2.x(:,1),temp.AlphaSweep{kk}.SOL_EN2.x(:,2),'color',[0.65,0.65,0.65])
    end
    subplot(1,2,1)
    plot(temp.SOL_LF2.x(:,1),temp.SOL_LF2.x(:,2),'LineWidth',2,'Color',[0.8,0,0.8])
    axis equal
    alpha_str{1} = num2str(temp.true_estgain/-20); 
    alpha_str{2} = 'Original';
    legend(alpha_str)
    xlim([-0.075, 0.025])
    xticks([-0.075:0.025:0.025])
    ylim([-0.14 0.14])
    yticks([-.1:0.05:.1])

    subplot(1,2,2)
    plot(temp.SOL_EN2.x(:,1),temp.SOL_EN2.x(:,2),'LineWidth',2,'Color',[0.8,0,0.8])
    axis equal
    alpha_str{1} = num2str(temp.true_estgain/-20); 
    alpha_str{2} = 'Original';
    legend(alpha_str)
    xlim([-0.025, 0.075])
    xticks([-0.025:0.025:0.075])
    ylim([-0.14 0.14])
    yticks([-.1:0.05:.1])
end
