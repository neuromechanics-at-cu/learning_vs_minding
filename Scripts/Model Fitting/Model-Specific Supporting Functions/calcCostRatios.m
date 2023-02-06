function [J,err2eff] = calcCostRatios(Params)

% Generate Trajectory
SOL = generateReachTrajectory_CurlOrBaseline(Params);

% Tracking Costs
J_LF2_tr_POSx = sum((SOL.x(1:end-1,1) - SOL.x(1:end-1,7)).^2*Params.Q(1,1));
J_LF2_tr_POSy = sum((SOL.x(1:end-1,2) - SOL.x(1:end-1,8)).^2*Params.Q(2,2));
J_LF2_tr_VEL = sum((SOL.x(1:end-1,3)).^2*Params.Q(3,3)) + ...
            sum((SOL.x(1:end-1,4)).^2*Params.Q(4,4));        
J_LF2_tr_F = sum((SOL.x(1:end-1,5)).^2*Params.Q(5,5)) + ...
            sum((SOL.x(1:end-1,6)).^2*Params.Q(6,6));    
J_LF2_tr_Fdot = sum((SOL.x(1:end-1,9)).^2*Params.Q(9,9)) + ...
            sum((SOL.x(1:end-1,10)).^2*Params.Q(10,10));     
J_LF2_tr_tot = J_LF2_tr_POSx + J_LF2_tr_POSy + J_LF2_tr_VEL + J_LF2_tr_F + J_LF2_tr_Fdot;

% End Costs
J_LF2_end_POS = (SOL.x(end,1) - SOL.x(end,7)).^2*Params.Phi(1,1) + ...
                (SOL.x(end,2) - SOL.x(end,8)).^2*Params.Phi(2,2);
J_LF2_end_VEL =  (SOL.x(end,3)).^2*Params.Phi(3,3) + ...
                (SOL.x(end,4)).^2*Params.Phi(4,4);
J_LF2_end_F =  (SOL.x(end,5)).^2*Params.Phi(5,5) + ...
                (SOL.x(end,6)).^2*Params.Phi(6,6);   
J_LF2_end_Fdot =  (SOL.x(end,9)).^2*Params.Phi(9,9) + ...
                (SOL.x(end,10)).^2*Params.Phi(10,10);  
J_LF2_end_tot = J_LF2_end_POS + J_LF2_end_VEL + J_LF2_end_F + J_LF2_end_Fdot;
% Control Costs
J_LF2_u = sum((SOL.u(1:end,1)).^2*Params.R(1,1)) + ...
            sum((SOL.u(1:end,2)).^2*Params.R(2,2)); 
% Total Cost
J_LF2_tot = J_LF2_tr_tot + J_LF2_end_tot + J_LF2_u;
% % Check Total Cost 
% J_LF2_tot_check = 0;
% for xxx = 1:size(SOL.x,1)-1
%     J_LF2_tot_check = J_LF2_tot_check + SOL.x(xxx,:)*Params.Q*SOL.x(xxx,:)' + SOL.u(xxx,:)*Params.R*SOL.u(xxx,:)';
% end
% J_LF2_tot_check = J_LF2_tot_check + SOL.x(end,:)*Params.Phi*SOL.x(end,:)';


% Find Total Costs, Proprotional Costs
% Tracking Costs
J_LF2_tr_POSx_sqState = (sum(SOL.x(1:end-1,1) - SOL.x(1:end-1,7)).^2)./(Params.trialLength-1);
J_LF2_tr_POSy_sqState = (sum(SOL.x(1:end-1,2) - SOL.x(1:end-1,8)).^2)./(Params.trialLength-1);
J_LF2_tr_VEL_sqState = (sum(SOL.x(1:end-1,3).^2) + ...
            sum(SOL.x(1:end-1,4).^2))./(Params.trialLength-1);        
J_LF2_tr_F_sqState = (sum(SOL.x(1:end-1,5).^2) + ...
            sum(SOL.x(1:end-1,6).^2))./(Params.trialLength-1);    
J_LF2_tr_Fdot_sqState = (sum(SOL.x(1:end-1,9).^2) + ...
            sum(SOL.x(1:end-1,10).^2))./(Params.trialLength-1);     
J_LF2_tr_tot_sqState = J_LF2_tr_POSx_sqState + J_LF2_tr_POSy_sqState + ...
                        J_LF2_tr_VEL_sqState + J_LF2_tr_F_sqState + J_LF2_tr_Fdot_sqState;
% End Costs
J_LF2_end_POS_sqState = (SOL.x(end,1) - SOL.x(end,7)).^2 + ...
                (SOL.x(end,2) - SOL.x(end,8)).^2;
J_LF2_end_VEL_sqState =  (SOL.x(end,3)).^2 + ...
                (SOL.x(end,4)).^2;
J_LF2_end_F_sqState =  (SOL.x(end,5)).^2 + ...
                (SOL.x(end,6)).^2;   
J_LF2_end_Fdot_sqState =  (SOL.x(end,9)).^2 + ...
                (SOL.x(end,10)).^2;  
J_LF2_end_tot_sqState = J_LF2_end_POS_sqState + J_LF2_end_VEL_sqState + J_LF2_end_F_sqState + J_LF2_end_Fdot_sqState;
% Control Costs
J_LF2_u_sqState = (sum(SOL.u(1:end,1).^2) + ...
            sum(SOL.u(1:end,2).^2))./(Params.trialLength-1); 
% Total Cost
J_LF2_tot_sqState = J_LF2_tr_tot_sqState*(Params.trialLength-1) + J_LF2_end_tot_sqState + J_LF2_u_sqState*(Params.trialLength-1);
% Check Total
J_LF2_tot_check = 0;
for xxx = 1:size(SOL.x,1)-1
    J_LF2_tot_check = J_LF2_tot_check + SOL.x(xxx,:)*SOL.x(xxx,:)' +  SOL.u(xxx,:)*SOL.u(xxx,:)';
end
J_LF2_tot_check = J_LF2_tot_check + SOL.x(end,:)*SOL.x(end,:)';


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


%% Just Cost Terms
% J_LF2_tr_POSx_norm = J_LF2_tr_POSx;
% J_LF2_tr_POSy_norm = J_LF2_tr_POSy;
% J_LF2_tr_VEL_norm = J_LF2_tr_VEL;
% J_LF2_tr_F_norm = J_LF2_tr_F;
% J_LF2_tr_Fdot_norm = J_LF2_tr_Fdot;
% J_LF2_end_POS_norm = J_LF2_end_POS;
% J_LF2_end_VEL_norm = J_LF2_end_VEL;
% J_LF2_end_F_norm = J_LF2_end_F;
% J_LF2_end_Fdot_norm = J_LF2_end_Fdot;
% J_LF2_u_norm = J_LF2_u;

%% Just Cost Weights
% J_LF2_tr_POSx_norm = Params.Q(1,1);
% J_LF2_tr_POSy_norm = Params.Q(2,2);
% J_LF2_tr_VEL_norm = Params.Q(3,3);
% J_LF2_tr_F_norm = Params.Q(5,5);
% J_LF2_tr_Fdot_norm = Params.Q(9,9);
% J_LF2_end_POS_norm = Params.Phi(1,1);
% J_LF2_end_VEL_norm = Params.Phi(3,3);
% J_LF2_end_F_norm = Params.Phi(5,5);
% J_LF2_end_Fdot_norm = Params.Phi(9,9);
% J_LF2_u_norm = Params.R(1,1);

% Combine Kinematic and Effort terms
% J_kin_combined = J_LF2_tr_POSx_norm+J_LF2_tr_POSy_norm+ J_LF2_tr_VEL_norm + J_LF2_end_POS_norm + J_LF2_end_VEL_norm;
% Removing velocity terms makes almost no difference, but helps disambiguate effort vs. error terms.

J_kin_combined = J_LF2_tr_POSx_norm+J_LF2_tr_POSy_norm+J_LF2_end_POS_norm; 
J_eff_combined = J_LF2_tr_F_norm+J_LF2_tr_Fdot_norm+J_LF2_u_norm + J_LF2_end_F_norm + J_LF2_end_Fdot_norm;

err2eff = J_kin_combined./J_eff_combined;

J.kin.J_LF2_tr_POSx_norm = J_LF2_tr_POSx_norm;
J.kin.J_LF2_tr_POSy_norm = J_LF2_tr_POSy_norm;
J.kin.J_LF2_end_POS_norm = J_LF2_end_POS_norm;

J.eff.J_LF2_tr_F_norm = J_LF2_tr_F_norm;
J.eff.J_LF2_tr_Fdot_norm = J_LF2_tr_Fdot_norm;
J.eff.J_LF2_u_norm = J_LF2_u_norm;
J.eff.J_LF2_end_F_norm = J_LF2_end_F_norm;
J.eff.J_LF2_end_Fdot_norm = J_LF2_end_Fdot_norm;
