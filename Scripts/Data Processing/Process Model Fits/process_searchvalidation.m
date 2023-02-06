%% Process Search Function Validation Fits

file_dir = '../../../Data/Model Fits/Reach Model and Search Validation (Fig 3)/All Fits';

CR_strings = {'HighCR','LowCR'};
learn_strings = {'0.6','1.0'};
searchtype_strings = {'fixed',... % Fixed Learning, Look for Costs
                      'fixedcostfindlearning',...
                      'varyall'};
numtraj_strings = {'singletraj','twotraj'};

for ii = 1:length(CR_strings)
    for jj = 1:length(learn_strings)
        for kk = 1:length(searchtype_strings)
            for ll = 1:length(numtraj_strings)
                load([file_dir,'/ModelValidation-',CR_strings{ii},'-',learn_strings{jj},searchtype_strings{kk},...
                    '-',numtraj_strings{ll},'.mat'],'AlphaSweep','Params','cost_ratios','SOL_LF2','SOL_EN2')
                % Add to Table
                % Single/Two Traj, Vary Costs?, Vary Learning?, True Cost Ratio, True Learning, ...
                % Found Mean Cost Ratio, Found STD Cost Ratio, Found Mean Learning, Found STD Learning
                
                % Save Relevant Parameters
                All_Searches{ii,jj,kk,ll}.AlphaSweep = AlphaSweep;
                All_Searches{ii,jj,kk,ll}.one_or_two_traj = ll;
                
                tParams = Params; tParams.estGain = -15; tParams.curlGain = 0;
                [~,true_costratio] = calcCostRatios(tParams);
                clear tParams
                All_Searches{ii,jj,kk,ll}.true_costratio = log(true_costratio);

                All_Searches{ii,jj,kk,ll}.true_estgain = Params.estGain;

                % If searching for proportion learned, save it
                if kk > 1
                    All_Searches{ii,jj,kk,ll}.vary_learning = true;
                    for nn = 1:10
                        estgains(nn) = AlphaSweep{nn}.Params.estGain;
                    end
                    All_Searches{ii,jj,kk,ll}.estgain_mean = mean(estgains);
                    All_Searches{ii,jj,kk,ll}.estgain_std = std(estgains);
                else
                    All_Searches{ii,jj,kk,ll}.vary_learning = false;
                    All_Searches{ii,jj,kk,ll}.estgain_mean = NaN;
                    All_Searches{ii,jj,kk,ll}.estgain_std = NaN;
                end
                % if searching over costs, save cost ratio
                if mod(kk,2) == 1
                    All_Searches{ii,jj,kk,ll}.vary_costs = true;
                    % Calculate using Learning / LF2
%                     All_Searches{ii,jj,kk,ll}.costratio_mean = mean(log(cost_ratios));
%                     All_Searches{ii,jj,kk,ll}.costratio_std = std(log(cost_ratios));

                    % Calculate using Washout/EN2 (A more accuate measure)
                    for nn = 1:10
                        tParams = AlphaSweep{nn}.Params;
                        tParams.curlGain = 0;
                        [~,newcost_ratios(nn)] = calcCostRatios(tParams);
                        clear tParams
                    end
                    All_Searches{ii,jj,kk,ll}.costratio_mean = mean(log(newcost_ratios));
                    All_Searches{ii,jj,kk,ll}.costratio_std = std(log(newcost_ratios));
                else
                    All_Searches{ii,jj,kk,ll}.vary_costs = false;
                    All_Searches{ii,jj,kk,ll}.costratio_mean = NaN;
                    All_Searches{ii,jj,kk,ll}.costratio_std = NaN;
                end
                % Save True Trajectories
                if exist('SOL_LF2','var')
                    All_Searches{ii,jj,kk,ll}.SOL_LF2 = SOL_LF2;
                end
                if exist('SOL_EN2','var')
                    All_Searches{ii,jj,kk,ll}.SOL_EN2 = SOL_EN2;
                end
            end
        end
    end
end

clear cost_ratios estgains Params AlphaSweep ii jj kk ll nn true_costratio ans newcost_ratios
save('../../../Data/Model Fits/Reach Model and Search Validation (Fig 3)/SearchFunctionValidationResults.mat')