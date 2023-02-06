% Model Validation Exercise 
% Produce Model Generated Trajectories
% Then Fit with a fixed alpha level

%% Generate Model Trajectory in a Curl Field
% Define Dynamics
x0 = [0;-0.1;0;0;0;0;0;0.1;0;0];
curlGain = -20;%Ns/m (curl gain)
dt = 1/200;%sec
Mass = [5 0;0 5];

invMass = inv(Mass);
t0 = 0;
Asys = [0 0 1 0 0 0 0 0 0 0; 
        0 0 0 1 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 1 0;
        0 0 0 0 0 0 0 0 0 1;
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 0 0 0];
Asys(3:4,5:6) = invMass;
Fsys = zeros(10);
Fsys(3:4,3:4) = invMass*[0 1;-1 0];
% Channel Trial Dynamics
stiffConst = 2000;
dampConst = 50;
Fsys_Chan = zeros(10);
Fsys_Chan(3:4,1:4) = invMass*[-1*stiffConst 0 -1*dampConst 0;
                               0            0  0           0];
% Controller Matrix
Bsys = zeros(10,2);
motorGain = 1e6;% to make R weights on same scale as Q and Phi
Bsys(9:10,1:2) = eye(2)*motorGain;
% Make discrete Time
AdisNull = eye(10) + Asys*dt;
AdisCurl = eye(10) + Asys*dt + Fsys*dt;
Bdis = dt*Bsys;

% Make Parameter Structure
Params.x0 = x0;
Params.dt = dt;
Params.Asys = Asys;
Params.B = Bdis;
Params.Fsys = Fsys;
Params.Fsys_Chan = Fsys_Chan;

%     Params.curlGain = curlGain;
Params.trialLength = 0.6*200;
Params.Mass = Mass;

Params.Q = zeros(size(Asys));
Params.R = eye(size(Bdis,2));
Params.Phi = eye(size(Asys));

Params.Mass = Mass;

% Plot Params
colors = [80,41,92;
          197,73,98;
          255,159,59]./255;
phase_names = {'Learning','Washout'};

%% For Learning and Washout
enviros = [-20,0];
prop_change = [0.2,1,5]; % ORIGINAL
alphas = 0:0.05:1;
        
% costRatios = nan(length(alphas),1);

for nn = 1:length(enviros)
    Params.curlGain = enviros(nn);
    
    for jj = 1:length(prop_change)
    
        q11 = 100*prop_change(jj); % ORIGINAL - FOR REFERENCE
        q22 = 10; % ORIGINAL
        q3344 = 0; % ORIGINAL
        q5566 = 1e-5; % ORIGINAL
        q991010 = 1e-10;
        r1122 = 5e2; % ORIGINAL - FOR REFERENCE
        phi1122 = 1e4*prop_change(jj); % ORIGINAL - FOR REFERENCE
        phi3344 = 1e-6;
        phi5566 = 1e-2;
        phi991010 = 1e-2;
        
        Params.Q(1,1) = q11; Params.Q(7,7) = q11;
        Params.Q(1,7) = -1*q11; Params.Q(7,1) = -1*q11;
        Params.Q(2,2) = q22; Params.Q(8,8) = q22;
        Params.Q(2,8) = -1*q22; Params.Q(8,2) = -1*q22;
        Params.Q(3,3) = q3344;
        Params.Q(5,5) = q5566; Params.Q(6,6) = q5566;
        Params.Q(9,9) = q991010; Params.Q(10,10) = q991010;
        
        Params.R(1,1) = r1122; Params.R(2,2) = r1122;
        
        Params.Phi(1,1) = phi1122; Params.Phi(7,7) = phi1122;
        Params.Phi(7,1) = -1*phi1122; Params.Phi(1,7) = -1*phi1122;
        Params.Phi(2,2) = phi1122; Params.Phi(8,8) = phi1122;
        Params.Phi(2,8) = -1*phi1122; Params.Phi(8,2) = -1*phi1122;
        Params.Phi(3,3) = phi3344; Params.Phi(4,4) = phi3344;
        Params.Phi(5,5) = phi5566; Params.Phi(6,6) = phi5566;
        Params.Phi(9,9) = phi991010; Params.Phi(10,10) = phi991010;
        
        %True CR
        Params.estGain = -15;
        [~,true_costratio(jj)] = calcCostRatios(Params);

        %% Plot these trajectories with different amounts of learning    
        for ii = 1:length(alphas)
            Params.estGain = alphas(ii)*curlGain;
            SOL = generateReachTrajectory_CurlOrBaseline(Params);
    
            figure(301+6*(nn-1)); hold on
            plot(SOL.x(:,1),SOL.x(:,2),'color',colors(jj,:),'LineWidth',2)
            [J_temp,crtemp] = calcCostRatios(Params);
            costRatios(ii,jj,nn) = crtemp; 
            J(ii,jj,nn) = J_temp;
    
            Metrics(ii,jj,nn) = calcErrorMetrics(Params);
            maxpx(ii,jj,nn) = Metrics(ii,jj,nn).MaxPx;
            maxfx(ii,jj,nn) = Metrics(ii,jj,nn).MaxFx*-1;
            adaptind(ii,jj,nn) = Metrics(ii,jj,nn).adaptInd/-20;

            SavedParams{ii,jj,nn} = Params;

            figure(302+6*(nn-1)); hold on
            plot(alphas(ii),Metrics(ii,jj,nn).MaxPx,'.','color',colors(jj,:),'MarkerSize',10)
            figure(303+6*(nn-1)); hold on
            plot(alphas(ii),Metrics(ii,jj,nn).MaxFx*-1,'.','color',colors(jj,:),'MarkerSize',10)
            figure(304+6*(nn-1)); hold on
            plot(alphas(ii),Metrics(ii,jj,nn).adaptInd/-20,'.','color',colors(jj,:),'MarkerSize',10)
        end
        figure(301+6*(nn-1))
        axis equal
        title(['Traj - ',phase_names{nn}])
        figure(302+6*(nn-1))
        title(['Max Px - ',phase_names{nn}])
        plot(alphas,maxpx(:,jj,nn),'color',colors(jj,:),'LineWidth',2)
        figure(303+6*(nn-1))
        title(['Max Fx - ',phase_names{nn}])
        plot(alphas,maxfx(:,jj,nn),'color',colors(jj,:),'LineWidth',2)
        ylim([0 20])
        figure(304+6*(nn-1))
        title(['Adapt Ind - ',phase_names{nn}])
        plot(alphas,adaptind(:,jj,nn),'color',colors(jj,:),'LineWidth',2)
        ylim([0 1])
        
        % Plot Cost Ratios and Each Cost Parameter
        if 0
            figure(305+6*(nn-1))
            hold on
            plot(alphas,log(costRatios(:,jj,nn)),'color',colors(jj,:),'LineWidth',2)
            plot(alphas,log(costRatios(:,jj,nn)),'.','color',colors(jj,:),'LineWidth',2,'MarkerSize',10)
            title(['Cost Ratio - ',phase_names{nn}])

            for ii = 1:length(alphas)
                figure(306+6*(nn-1))
                subplot(2,4,1)
                hold on
                plot(alphas(ii),J(ii,jj,nn).kin.J_LF2_tr_POSx_norm,'x','color',colors(jj,:))
                title('Kin: x pos')
            
                subplot(2,4,2)
                hold on
                plot(alphas(ii),J(ii,jj,nn).kin.J_LF2_tr_POSy_norm,'x','color',colors(jj,:))
                title('Kin: y pos')
            
                subplot(2,4,3)
                hold on
                plot(alphas(ii),J(ii,jj,nn).kin.J_LF2_end_POS_norm,'x','color',colors(jj,:))
                title('Kin: End Pos')
            
                subplot(2,4,4)
                hold on
                plot(alphas(ii),J(ii,jj,nn).eff.J_LF2_tr_F_norm,'x','color',colors(jj,:))
                title('Eff: F')
            
                subplot(2,4,5)
                hold on
                plot(alphas(ii),J(ii,jj,nn).eff.J_LF2_tr_Fdot_norm,'x','color',colors(jj,:))
                title('Eff: Fdot')
            
                subplot(2,4,6)
                hold on
                plot(alphas(ii),J(ii,jj,nn).eff.J_LF2_u_norm,'x','color',colors(jj,:))
                title('Eff: U/Fddot')
            
                subplot(2,4,7)
                hold on
                plot(alphas(ii),J(ii,jj,nn).eff.J_LF2_end_F_norm,'x','color',colors(jj,:))
                title('Eff: end F')
            
                subplot(2,4,8)
                hold on
                plot(alphas(ii),J(ii,jj,nn).eff.J_LF2_end_Fdot_norm,'x','color',colors(jj,:))
                title('Eff: end Fdot')
            end
            sgtitle(phase_names{nn})
        end
    end
end

%% Match Simple Model Predictions
alpha_index_to_plot = [17,15,13,11];

types_of_lines = {'-','--',':','-.'};

figure(313)
hold on
for ii = 1:length(alpha_index_to_plot) 
    true_px(ii) = maxpx(alpha_index_to_plot(ii),2,1);
    learn_cr1(ii) = interp1(maxpx(:,1,1),alphas,true_px(ii),'linear');
    learn_cr3(ii) = interp1(maxpx(:,3,1),alphas,true_px(ii),'linear');
    
    plot([learn_cr1(ii),alphas(alpha_index_to_plot(ii)),learn_cr3(ii)],log(costRatios(1,:,1)),'k','LineWidth',2,'LineStyle',types_of_lines{ii})
    plot(learn_cr1(ii),log(costRatios(1,1,1)),'.','Color',colors(1,:),'MarkerSize',20)
    plot(alphas(alpha_index_to_plot(ii)),log(costRatios(1,2,1)),'.','Color',colors(2,:),'MarkerSize',20)
    plot(learn_cr3(ii),log(costRatios(1,3,1)),'.','Color',colors(3,:),'MarkerSize',20)
end
xlim([0.2 1])
title('Iso-Maximum Horizontal Error Lines')

% % Plot iso-maxpx Trajectories 
% for ii = 1:length(alpha_index_to_plot) 
%     figure(319+ii)
%     hold on
%     SavedParams{1,1}.estGain = learn_cr1(ii)*-20;
%     SOL_cr1 = generateReachTrajectory_CurlOrBaseline(SavedParams{1,1});
%     plot(SOL_cr1.x(:,1),SOL_cr1.x(:,2),'color',colors(1,:))
% 
%     SavedParams{1,2}.estGain = alphas(alpha_index_to_plot(ii))*-20;
%     SOL_cr2 = generateReachTrajectory_CurlOrBaseline(SavedParams{1,2});
%     plot(SOL_cr2.x(:,1),SOL_cr2.x(:,2),'color',colors(2,:))
% 
%     SavedParams{1,3}.estGain = learn_cr3(ii)*-20;
%     SOL_cr3 = generateReachTrajectory_CurlOrBaseline(SavedParams{1,3});
%     plot(SOL_cr3.x(:,1),SOL_cr3.x(:,2),'color',colors(3,:))
%     
%     axis equal
% end
   
%% Plot Other Learning Metrics
if 0
    %% Adaptation Index
    figure(314)
    hold on
    for ii = 1:length(alpha_index_to_plot) 
        true_adaptind(ii) = adaptind(alpha_index_to_plot(ii),2,1);
        learn_cr1(ii) = interp1(adaptind(:,1,1),alphas,true_adaptind(ii),'linear');
        learn_cr3(ii) = interp1(adaptind(:,3,1),alphas,true_adaptind(ii),'linear');
        
        plot([learn_cr1(ii),alphas(alpha_index_to_plot(ii)),learn_cr3(ii)],log(costRatios(1,:,1)),'k','LineWidth',2,'LineStyle',types_of_lines{ii})
        plot(learn_cr1(ii),log(costRatios(1,1,1)),'.','Color',colors(1,:),'MarkerSize',20)
        plot(alphas(alpha_index_to_plot(ii)),log(costRatios(1,2,1)),'.','Color',colors(2,:),'MarkerSize',20)
        plot(learn_cr3(ii),log(costRatios(1,3,1)),'.','Color',colors(3,:),'MarkerSize',20)
    end
    xlim([0.2 1])
    title('Iso-Adaptation Index Lines')
    
    
    % % Plot Trajectories 
    % for ii = 1:length(alpha_index_to_plot) 
    %     figure(339+ii)
    %     hold on
    %     SavedParams{1,1}.estGain = learn_cr1(ii)*-20;
    %     SOL_cr1 = generateReachTrajectory_CurlOrBaseline(SavedParams{1,1});
    %     plot(SOL_cr1.x(:,1),SOL_cr1.x(:,2),'color',colors(1,:))
    % 
    %     SavedParams{1,2}.estGain = alphas(alpha_index_to_plot(ii))*-20;
    %     SOL_cr2 = generateReachTrajectory_CurlOrBaseline(SavedParams{1,2});
    %     plot(SOL_cr2.x(:,1),SOL_cr2.x(:,2),'color',colors(2,:))
    % 
    %     SavedParams{1,3}.estGain = learn_cr3(ii)*-20;
    %     SOL_cr3 = generateReachTrajectory_CurlOrBaseline(SavedParams{1,3});
    %     plot(SOL_cr3.x(:,1),SOL_cr3.x(:,2),'color',colors(3,:))
    %     
    %     axis equal
    % end
    
    %% Max Fx
    figure(315)
    hold on
    for ii = 1:length(alpha_index_to_plot) 
        true_maxfx(ii) = maxfx(alpha_index_to_plot(ii),2,1);
        learn_cr1(ii) = interp1(maxfx(:,1,1),alphas,true_maxfx(ii),'linear');
        learn_cr3(ii) = interp1(maxfx(:,3,1),alphas,true_maxfx(ii),'linear');
        
        plot([learn_cr1(ii),alphas(alpha_index_to_plot(ii)),learn_cr3(ii)],log(costRatios(1,:,1)),'k','LineWidth',2,'LineStyle',types_of_lines{ii})
        plot(learn_cr1(ii),log(costRatios(1,1,1)),'.','Color',colors(1,:),'MarkerSize',20)
        plot(alphas(alpha_index_to_plot(ii)),log(costRatios(1,2,1)),'.','Color',colors(2,:),'MarkerSize',20)
        plot(learn_cr3(ii),log(costRatios(1,3,1)),'.','Color',colors(3,:),'MarkerSize',20)
    end
    xlim([0.2 1])
    title('Iso-Maximum Horizontal Force Lines')
    
    % % Plot Trajectories 
    % for ii = 1:length(alpha_index_to_plot) 
    %     figure(359+ii)
    %     hold on
    %     SavedParams{1,1}.estGain = learn_cr1(ii)*-20;
    %     SOL_cr1 = generateReachTrajectory_CurlOrBaseline(SavedParams{1,1});
    %     plot(SOL_cr1.x(:,1),SOL_cr1.x(:,2),'color',colors(1,:))
    % 
    %     SavedParams{1,2}.estGain = alphas(alpha_index_to_plot(ii))*-20;
    %     SOL_cr2 = generateReachTrajectory_CurlOrBaseline(SavedParams{1,2});
    %     plot(SOL_cr2.x(:,1),SOL_cr2.x(:,2),'color',colors(2,:))
    % 
    %     SavedParams{1,3}.estGain = learn_cr3(ii)*-20;
    %     SOL_cr3 = generateReachTrajectory_CurlOrBaseline(SavedParams{1,3});
    %     plot(SOL_cr3.x(:,1),SOL_cr3.x(:,2),'color',colors(3,:))
    %     
    %     axis equal
    % end
end