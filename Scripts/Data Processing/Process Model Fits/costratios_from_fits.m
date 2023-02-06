%% Calculate Costs from Model Fits

% Data to Load
subjNames = {'Y','O'};
runNames = 1:10;
% runNames = 1:9;

% alphaNames = 0.6:0.05:0.8;
alphaNames = 0.65:0.05:0.85;

% colors = [1,0,0;1,127/255,0];
colors = [0,128/255,0;0,0,128/255];

% curDir = pwd;

% Initialize Cost Matrices
Q_11 = nan(length(subjNames),length(runNames),length(alphaNames));
Q_22 = nan(length(subjNames),length(runNames),length(alphaNames));
Q_3344 = nan(length(subjNames),length(runNames),length(alphaNames));
Q_5566 = nan(length(subjNames),length(runNames),length(alphaNames));
Q_991010 = nan(length(subjNames),length(runNames),length(alphaNames));
R_1122 = nan(length(subjNames),length(runNames),length(alphaNames));
PHI_1122 = nan(length(subjNames),length(runNames),length(alphaNames));
PHI_3344 = nan(length(subjNames),length(runNames),length(alphaNames));
PHI_5566 = nan(length(subjNames),length(runNames),length(alphaNames));
PHI_991010 = nan(length(subjNames),length(runNames),length(alphaNames));

% Calculate total costs by state
for iii = 1:length(subjNames)
    for jjj = 1:length(alphaNames)
        for kkk = 1:length(runNames)
            if exist(['../../../Data/Model Fits/Cost Weight Analysis Fits (Fig 6)/All Fits/',...
                    subjNames{iii},'-',num2str(runNames(kkk)),'-alpha',num2str(alphaNames(jjj)),'-1.96.mat']) %TODO: generalized to CI
                load(['../../../Data/Model Fits/Cost Weight Analysis Fits (Fig 6)/All Fits/',...
                    subjNames{iii},'-',num2str(runNames(kkk)),'-alpha',num2str(alphaNames(jjj)),'-1.96.mat'])  %TODO: generalized to CI
                trialLength(iii) = length(yData.LF2.PxSubjMean);
                Q_11(iii,kkk,jjj) = Young.LF2.Params.Q(1,1);
                Q_22(iii,kkk,jjj) = Young.LF2.Params.Q(2,2);
                Q_3344(iii,kkk,jjj) = Young.LF2.Params.Q(3,3);
                Q_5566(iii,kkk,jjj) = Young.LF2.Params.Q(5,5);
                Q_991010(iii,kkk,jjj) = Young.LF2.Params.Q(9,9);
                R_1122(iii,kkk,jjj) = Young.LF2.Params.R(1,1);
                PHI_1122(iii,kkk,jjj) = Young.LF2.Params.Phi(1,1);
                PHI_3344(iii,kkk,jjj) = Young.LF2.Params.Phi(3,3);
                PHI_5566(iii,kkk,jjj) = Young.LF2.Params.Phi(5,5);
                PHI_991010(iii,kkk,jjj) = Young.LF2.Params.Phi(9,9);
                trialLength(iii) = length(yData.LF2.PxSubjMean);
  
            % Find Total Costs, Proprotional Costs
                % Tracking Costs
                J_LF2_tr_POSx(iii,kkk,jjj) = sum((SOL_LF2.x(1:end-1,1) - SOL_LF2.x(1:end-1,7)).^2*Young.LF2.Params.Q(1,1));
                J_LF2_tr_POSy(iii,kkk,jjj) = sum((SOL_LF2.x(1:end-1,2) - SOL_LF2.x(1:end-1,8)).^2*Young.LF2.Params.Q(2,2));
                J_LF2_tr_VEL(iii,kkk,jjj) = sum((SOL_LF2.x(1:end-1,3)).^2*Young.LF2.Params.Q(3,3)) + ...
                                sum((SOL_LF2.x(1:end-1,4)).^2*Young.LF2.Params.Q(4,4));        
                J_LF2_tr_F(iii,kkk,jjj) = sum((SOL_LF2.x(1:end-1,5)).^2*Young.LF2.Params.Q(5,5)) + ...
                                sum((SOL_LF2.x(1:end-1,6)).^2*Young.LF2.Params.Q(6,6));    
                J_LF2_tr_Fdot(iii,kkk,jjj) = sum((SOL_LF2.x(1:end-1,9)).^2*Young.LF2.Params.Q(9,9)) + ...
                                sum((SOL_LF2.x(1:end-1,10)).^2*Young.LF2.Params.Q(10,10));     
                J_LF2_tr_tot(iii,kkk,jjj) = J_LF2_tr_POSx(iii,kkk,jjj) + J_LF2_tr_POSy(iii,kkk,jjj) + ...
                                            J_LF2_tr_VEL(iii,kkk,jjj) + J_LF2_tr_F(iii,kkk,jjj) + J_LF2_tr_Fdot(iii,kkk,jjj);

                % End Costs
                J_LF2_end_POS(iii,kkk,jjj) = (SOL_LF2.x(end,1) - SOL_LF2.x(end,7)).^2*Young.LF2.Params.Phi(1,1) + ...
                                    (SOL_LF2.x(end,2) - SOL_LF2.x(end,8)).^2*Young.LF2.Params.Phi(2,2);
                J_LF2_end_VEL(iii,kkk,jjj) =  (SOL_LF2.x(end,3)).^2*Young.LF2.Params.Phi(3,3) + ...
                                    (SOL_LF2.x(end,4)).^2*Young.LF2.Params.Phi(4,4);
                J_LF2_end_F(iii,kkk,jjj) =  (SOL_LF2.x(end,5)).^2*Young.LF2.Params.Phi(5,5) + ...
                                    (SOL_LF2.x(end,6)).^2*Young.LF2.Params.Phi(6,6);   
                J_LF2_end_Fdot(iii,kkk,jjj) =  (SOL_LF2.x(end,9)).^2*Young.LF2.Params.Phi(9,9) + ...
                                    (SOL_LF2.x(end,10)).^2*Young.LF2.Params.Phi(10,10);  
                J_LF2_end_tot(iii,kkk,jjj) = J_LF2_end_POS(iii,kkk,jjj) + J_LF2_end_VEL(iii,kkk,jjj) + J_LF2_end_F(iii,kkk,jjj) + J_LF2_end_Fdot(iii,kkk,jjj);
                % Control Costs
                J_LF2_u(iii,kkk,jjj) = sum((SOL_LF2.u(1:end,1)).^2*Young.LF2.Params.R(1,1)) + ...
                                sum((SOL_LF2.u(1:end,2)).^2*Young.LF2.Params.R(2,2)); 
                % Total Cost
                J_LF2_tot(iii,kkk,jjj) = J_LF2_tr_tot(iii,kkk,jjj) + J_LF2_end_tot(iii,kkk,jjj) + J_LF2_u(iii,kkk,jjj);
                % Check Total Cost 
                J_LF2_tot_check(iii,kkk,jjj) = 0;
                for xxx = 1:size(SOL_LF2.x,1)-1
                    J_LF2_tot_check(iii,kkk,jjj) = J_LF2_tot_check(iii,kkk,jjj) + SOL_LF2.x(xxx,:)*Young.LF2.Params.Q*SOL_LF2.x(xxx,:)' +  SOL_LF2.u(xxx,:)*Young.LF2.Params.R*SOL_LF2.u(xxx,:)';
                end
                J_LF2_tot_check(iii,kkk,jjj) = J_LF2_tot_check(iii,kkk,jjj) + SOL_LF2.x(end,:)*Young.LF2.Params.Phi*SOL_LF2.x(end,:)';

                % Find Total Costs, Proprotional Costs
                % Tracking Costs
                J_LF2_tr_POSx_sqState(iii,kkk,jjj) = (sum(SOL_LF2.x(1:end-1,1) - SOL_LF2.x(1:end-1,7)).^2)./(trialLength(iii)-1);
                J_LF2_tr_POSy_sqState(iii,kkk,jjj) = (sum(SOL_LF2.x(1:end-1,2) - SOL_LF2.x(1:end-1,8)).^2)./(trialLength(iii)-1);
                J_LF2_tr_VEL_sqState(iii,kkk,jjj) = (sum(SOL_LF2.x(1:end-1,3).^2) + ...
                                sum(SOL_LF2.x(1:end-1,4).^2))./(trialLength(iii)-1);        
                J_LF2_tr_F_sqState(iii,kkk,jjj) = (sum(SOL_LF2.x(1:end-1,5).^2) + ...
                                sum(SOL_LF2.x(1:end-1,6).^2))./(trialLength(iii)-1);    
                J_LF2_tr_Fdot_sqState(iii,kkk,jjj) = (sum(SOL_LF2.x(1:end-1,9).^2) + ...
                                sum(SOL_LF2.x(1:end-1,10).^2))./(trialLength(iii)-1);     
                J_LF2_tr_tot_sqState(iii,kkk,jjj) = J_LF2_tr_POSx_sqState(iii,kkk,jjj) + J_LF2_tr_POSy_sqState(iii,kkk,jjj) + ...
                                            J_LF2_tr_VEL_sqState(iii,kkk,jjj) + J_LF2_tr_F_sqState(iii,kkk,jjj) + J_LF2_tr_Fdot_sqState(iii,kkk,jjj);
                % End Costs
                J_LF2_end_POS_sqState(iii,kkk,jjj) = (SOL_LF2.x(end,1) - SOL_LF2.x(end,7)).^2 + ...
                                    (SOL_LF2.x(end,2) - SOL_LF2.x(end,8)).^2;
                J_LF2_end_VEL_sqState(iii,kkk,jjj) =  (SOL_LF2.x(end,3)).^2 + ...
                                    (SOL_LF2.x(end,4)).^2;
                J_LF2_end_F_sqState(iii,kkk,jjj) =  (SOL_LF2.x(end,5)).^2 + ...
                                    (SOL_LF2.x(end,6)).^2;   
                J_LF2_end_Fdot_sqState(iii,kkk,jjj) =  (SOL_LF2.x(end,9)).^2 + ...
                                    (SOL_LF2.x(end,10)).^2;  
                J_LF2_end_tot_sqState(iii,kkk,jjj) = J_LF2_end_POS_sqState(iii,kkk,jjj) + J_LF2_end_VEL_sqState(iii,kkk,jjj) + J_LF2_end_F_sqState(iii,kkk,jjj) + J_LF2_end_Fdot_sqState(iii,kkk,jjj);
                % Control Costs
                J_LF2_u_sqState(iii,kkk,jjj) = (sum(SOL_LF2.u(1:end,1).^2) + ...
                                sum(SOL_LF2.u(1:end,2).^2))./(trialLength(iii)-1); 
                % Total Cost
                J_LF2_tot_sqState(iii,kkk,jjj) = J_LF2_tr_tot_sqState(iii,kkk,jjj)*(trialLength(iii)-1) + J_LF2_end_tot_sqState(iii,kkk,jjj) + J_LF2_u_sqState(iii,kkk,jjj)*(trialLength(iii)-1);
                % Check Total
                J_LF2_tot_check(iii,kkk,jjj) = 0;
                for xxx = 1:size(SOL_LF2.x,1)-1
                    J_LF2_tot_check(iii,kkk,jjj) = J_LF2_tot_check(iii,kkk,jjj) + SOL_LF2.x(xxx,:)*SOL_LF2.x(xxx,:)' +  SOL_LF2.u(xxx,:)*SOL_LF2.u(xxx,:)';
                end
                J_LF2_tot_check(iii,kkk,jjj) = J_LF2_tot_check(iii,kkk,jjj) + SOL_LF2.x(end,:)*SOL_LF2.x(end,:)';
            else
               disp(['Cant find trial ',num2str(iii),' ',num2str(jjj),' ',num2str(kkk)])
            end
        end
    end
end

%% Normalize and Combine Cost Terms
J_LF2_tr_POSx_norm = J_LF2_tr_POSx./J_LF2_tr_POSx_sqState;
J_LF2_tr_POSy_norm = J_LF2_tr_POSy./J_LF2_tr_POSy_sqState;
J_LF2_tr_VEL_norm = J_LF2_tr_VEL./J_LF2_tr_VEL_sqState;
J_LF2_tr_F_norm = J_LF2_tr_F./J_LF2_tr_F_sqState;
J_LF2_tr_Fdot_norm = J_LF2_tr_Fdot./J_LF2_tr_Fdot_sqState;
J_LF2_end_POS_norm = J_LF2_end_POS./J_LF2_end_POS_sqState;
J_LF2_end_VEL_norm = J_LF2_end_VEL./J_LF2_end_VEL_sqState;
J_LF2_end_F_norm = J_LF2_end_F./J_LF2_end_F_sqState;
J_LF2_end_Fdot_norm = J_LF2_end_Fdot./J_LF2_end_Fdot_sqState;
J_LF2_u_norm = J_LF2_u./J_LF2_u_sqState;

% Combine Kinematic and Effort terms
% J_kin_combined = J_LF2_tr_POSx_norm+J_LF2_tr_POSy_norm+ J_LF2_tr_VEL_norm + J_LF2_end_POS_norm + J_LF2_end_VEL_norm;
% Removing velocity terms makes almost no difference, but helps disambiguate effort vs. error terms.
J_kin_combined = J_LF2_tr_POSx_norm+J_LF2_tr_POSy_norm+J_LF2_end_POS_norm; 
J_eff_combined = J_LF2_tr_F_norm+J_LF2_tr_Fdot_norm+J_LF2_u_norm + J_LF2_end_F_norm + J_LF2_end_Fdot_norm;

err2force = J_kin_combined./J_eff_combined;

%% Stats
% Reshape
% alphaMat = repmat(alphaNames,10,1);
% vecAlpha = reshape(alphaMat,1,50);
% vecErr2Force_Y = reshape(squeeze(err2force(1,:,:)),1,50);
% vecErr2Force_O = reshape(squeeze(err2force(2,:,:)),1,50);
alphaMat = repmat(alphaNames,length(runNames),1);
vecAlpha = reshape(alphaMat,1,length(alphaNames)*length(runNames));
vecErr2Force_Y = reshape(squeeze(err2force(1,:,:)),1,length(alphaNames)*length(runNames));
vecErr2Force_O = reshape(squeeze(err2force(2,:,:)),1,length(alphaNames)*length(runNames));

% [B_Y,BINT_Y,R_Y,RINT_Y,STATS_Y] = regress(vecErr2Force_Y',[vecAlpha',ones(size(vecErr2Force_Y'))]);
% [B_O,BINT_O,R_O,RINT_O,STATS_O] = regress(vecErr2Force_O',[vecAlpha',ones(size(vecErr2Force_O'))]);

[B_Y_log,BINT_Y_log,R_Y_log,RINT_Y_log,STATS_Y_log] = regress(log(vecErr2Force_Y)',[vecAlpha',ones(size(vecErr2Force_Y'))]);
[B_O_log,BINT_O_log,R_O_log,RINT_O_log,STATS_O_log] = regress(log(vecErr2Force_O)',[vecAlpha',ones(size(vecErr2Force_O'))]);

%% p-values for slope
p_Y = STATS_Y_log(3);
p_O = STATS_O_log(3);

%% Confidence Intervals of Samples
% Bootstrapping
ci_Y_log = bootci(10000,@mean,log(vecErr2Force_Y));
ci_O_log = bootci(10000,@mean,log(vecErr2Force_O));

save('../../../Data/Model Fits/Cost Weight Analysis Fits (Fig 6)/CostRatios_From_Fits.mat')