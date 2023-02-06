% Model Validation Continued
% Combine LF2 and EN2

%% Run Model Fitting Code

%% Generate Trajectory to be fit
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

Params.curlGain = curlGain;
Params.trialLength = 0.6*200;
Params.Mass = Mass;

Params.Q = zeros(size(Asys));
Params.R = eye(size(Bdis,2));
Params.Phi = eye(size(Asys));

Params.Mass = Mass;

%% Specify Costs
prop_cost = [5,1,0.2];

q11 = 100*prop_cost(cost_ratio_type); 
q22 = 10;
q3344 = 0;
q5566 = 1e-5;
q991010 = 0;
r1122 = 5e2;
phi1122 = 1e4*prop_cost(cost_ratio_type); 
phi3344 = 1e-6; 
phi5566 = 1e-2; 
phi991010 = 1e-2;

%Assign Base Values to Initial Params
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


%%%%% TRUE LEVEL OF LEARNING %%%%%%%
Params.estGain = prop_learned*curlGain;

%LF2
SOL_LF2 = generateReachTrajectory_CurlOrBaseline(Params);
SOL_LF2_Chan = generateReachTrajectory_Channel(Params);
stiffConst = 2000;
dampConst = 50;
ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_LF2_Chan.x(:,1)';SOL_LF2_Chan.x(:,2)';SOL_LF2_Chan.x(:,3)';SOL_LF2_Chan.x(:,4)'];
SOL_LF2_Chan.Fx = ForceMat(1,:)';
SOL_LF2_Chan.Fy = ForceMat(2,:)';
[~,mfx_ind] = max(abs(SOL_LF2_Chan.Fx));
SOL_LF2_Chan.FxMax = SOL_LF2_Chan.Fx(mfx_ind);
SOL_LF2_Chan.AdaptInd = fitgain(SOL_LF2_Chan.Fx,SOL_LF2_Chan.x(:,4)); 

%EN2
Params.curlGain = 0;
SOL_EN2 = generateReachTrajectory_CurlOrBaseline(Params);
SOL_EN2_Chan = generateReachTrajectory_Channel(Params);
stiffConst = 2000;
dampConst = 50;
ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_EN2_Chan.x(:,1)';SOL_EN2_Chan.x(:,2)';SOL_EN2_Chan.x(:,3)';SOL_EN2_Chan.x(:,4)'];
SOL_EN2_Chan.Fx = ForceMat(1,:)';
SOL_EN2_Chan.Fy = ForceMat(2,:)';
[~,mfx_ind] = max(abs(SOL_EN2_Chan.Fx));
SOL_EN2_Chan.FxMax = SOL_EN2_Chan.Fx(mfx_ind);
SOL_EN2_Chan.AdaptInd = fitgain(SOL_EN2_Chan.Fx,SOL_EN2_Chan.x(:,4)); 
Params.curlGain = -20;

%% Assign What Parameters to Vary for Model Fitting

% In this case, it is all cost weights, but the proportion learned is fixed
FitIndices = [vary_costs, ... % R(1,1) = R(2,2) (control effort) (second deriv of force)
              vary_costs, ... % Q(5,5) = Q(6,6) (force penalty)
              vary_costs, ... % Q(1,1) = Q(7,7) = -Q(1,7) = -Q(7,1) (horizontal tracking)             
              vary_costs, ... % Q(2,2) = Q(8,8) = -Q(2,8) = -Q(8,2) (vertical dist. penalty)             
              0, ... % Q(3,3) % Horizontal Velocity
              vary_costs, ... % Phi(1,1) = Phi(2,2) = Phi(7,7) = Phi(8,8) = -Phi(1,7)... (terminal position)
              vary_costs, ... % Phi(3,3) = Phi(4,4) (terminal velocity)
              vary_costs, ... % Phi(5,5) = Phi(6,6) (terminal force)
              vary_learning, ... % estGain
              vary_costs, ... % Q(9,9) = Q(10,10) (rate of force)
              vary_costs];    % Phi(9,9) = Phi(10,10) (terminal rate of force)

FitIndices = logical(FitIndices);
numX = sum(FitIndices);
numsols = 10;
AA = [];
bb = [];
Aeq = [];
beq = [];
lb = zeros(length(FitIndices),1);
lb(1) = 1;
lb(6) = 1e-20;
lb(9) = -20*1.4;
lb = lb(FitIndices);
ub = 1e20*ones(length(FitIndices),1);
ub(9) = 0;
ub = ub(FitIndices);
nonlcon = [];
options = optimset('MaxFunEval',1E5,'MaxIter',1E5,'TolFun',1E-14,...
        'TolX',1E-14,'TolCon',1E-14,'LargeScale','on','Algorithm','sqp');

Values = zeros(numsols,1);
Results = zeros(numsols,numX);
Flags = zeros(numsols,1);
% Initial FMINCON Parameters
%            [R(1,1), Q(5,5), Q(1,1), Q(2,2), Q(3,3), Phi(1,1), Phi(3,3), Phi(5,5),
Starting_X = [r1122;  q5566;  q11;    q22;    q3344;  phi1122;  phi3344;  phi5566;      ...
                %estGain, rate of Force; terminal Rate of Force]
                Params.estGain;   0;  0 ]; 


for jj = 1:length(alphas)
    bestFitParams = Params;
    bestFitParams.estGain = alphas(jj)*curlGain;
    kk = 1;
    while kk <= n_sols
        alphas(jj)
        kk
%         bestFitParams = Params;
        Values = zeros(numsols,1);
        Results = zeros(numsols,numX);
        Flags = zeros(numsols,1);

        
        for m = 1:numsols 
            disp(['minimization loop: ',num2str(m)])
    
            Xo = Starting_X;
            % Randomize Starting Point
            Xo(1:8) = Xo(1:8) + 10*rand(8,1).*Xo(1:8);
            Xo(9) = Xo(9) + 2*(rand(1,1)-0.5).*Xo(9);
            Xo(10:11) = Xo(10:11) + 10*rand(2,1).*Xo(10:11);
    
            Xo = Xo(FitIndices);
    
%             [X, fval, flag] = fmincon(@(X) qualityOfSingleTrajectory_noCI(X,SOL,SOL_Chan,FitIndices,bestFitParams),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
            [X, fval, flag] = fmincon(@(X) qualityOfLF2EN2Trajectory_noCI(X,SOL_LF2,SOL_LF2_Chan,SOL_EN2,SOL_EN2_Chan,FitIndices,bestFitParams),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
            Results(m,:) = X';
            Values(m) = fval;
            Flags(m) = flag;
        end
    
        % Refined Search
    %     options = optimset('MaxFunEval',1E4,'MaxIter',1E4,'TolFun',1E-18,...
    %         'TolX',1E-16,'TolCon',1E-16,'LargeScale','on','Algorithm','sqp');
        options = optimset('MaxFunEval',1E5,'MaxIter',1E5,'TolFun',1E-18,...
            'TolX',1E-16,'TolCon',1E-16,'LargeScale','on','Algorithm','sqp');
    
        for m = numsols+1:2*numsols
            disp(['refining minimization loop: ',num2str(m)])
            Xo = Results(m-numsols,:)';
%             [X, fval, flag] = fmincon(@(X) qualityOfSingleTrajectory_noCI(X,SOL,SOL_Chan,FitIndices,bestFitParams),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
            [X, fval, flag] = fmincon(@(X) qualityOfLF2EN2Trajectory_noCI(X,SOL_LF2,SOL_LF2_Chan,SOL_EN2,SOL_EN2_Chan,FitIndices,bestFitParams),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);

            Results(m,:) = X';
            Values(m) = fval;
            Flags(m) = flag;
        end
        
        % Add bounding fits to comparison
        Results(2*numsols+1,:) = lb';
        Results(2*numsols+1,end) = ub(end);
%         Values(2*numsols+1) = qualityOfSingleTrajectory_noCI(lb,SOL,SOL_Chan,FitIndices,bestFitParams);
        Values(2*numsols+1) = qualityOfLF2EN2Trajectory_noCI(lb,SOL_LF2,SOL_LF2_Chan,SOL_EN2,SOL_EN2_Chan,FitIndices,bestFitParams);
    
        [val, in] = min(Values);
        BestX = Results(in,:);
    
        if in == 2*numsols + 1
            pause
            [val,in] = min(Values(1:end-1));
            BestX = Results(in,:);
        end
    
        % Generate Trajectory from Solution
%         tParams = Params;
        tParams = bestFitParams;
        count = 1;
        if FitIndices(1)
            tParams.R(1,1) = BestX(count);
            tParams.R(2,2) = BestX(count);
            count = count + 1;
        end
        if FitIndices(2)
            tParams.Q(5,5) = BestX(count);
            tParams.Q(6,6) = BestX(count);
            count = count + 1;
        end
        if FitIndices(3)
            tParams.Q(1,1) = BestX(count);
            tParams.Q(7,7) = BestX(count);
            tParams.Q(1,7) = -1*BestX(count);
            tParams.Q(7,1) = -1*BestX(count);
            count = count + 1;
        end
        if FitIndices(4) 
            tParams.Q(2,2) = BestX(count);
            tParams.Q(8,8) = BestX(count);
            tParams.Q(2,8) = -1*BestX(count);
            tParams.Q(8,2) = -1*BestX(count);
            count = count + 1;
        end
        if FitIndices(5)
           tParams.Q(3,3) = BestX(count); 
           count = count + 1;
        end
        if FitIndices(6)
            tParams.Phi(1,1) = BestX(count);
            tParams.Phi(7,7) = BestX(count);
            tParams.Phi(1,7) = -1*BestX(count);
            tParams.Phi(7,1) = -1*BestX(count);
            tParams.Phi(2,2) = BestX(count);
            tParams.Phi(8,8) = BestX(count);
            tParams.Phi(2,8) = -1*BestX(count);
            tParams.Phi(8,2) = -1*BestX(count);
            count = count + 1;
        end
        if FitIndices(7)
            tParams.Phi(3,3) = BestX(count);
            tParams.Phi(4,4) = BestX(count);
            count = count + 1;
        end
        if FitIndices(8)
            tParams.Phi(5,5) = BestX(count);
            tParams.Phi(6,6) = BestX(count);
            count = count + 1;
        end
        if FitIndices(9)
            tParams.estGain = BestX(count);
            count = count + 1;
%         else
%             tParams.estGain = curlGain*alphas(jj);
        end
        if FitIndices(10)
            tParams.Q(9,9) = BestX(count);
            tParams.Q(10,10) = BestX(count);
            count = count + 1;
        end
        if FitIndices(11)
            tParams.Phi(9,9) = BestX(count);
            tParams.Phi(10,10) = BestX(count);
        end
%         SOL_new = generateReachTrajectory_CurlOrBaseline(tParams);
        SOL_LF2_new = generateReachTrajectory_CurlOrBaseline(tParams);
        tParams.curlGain = 0;
        SOL_EN2_new = generateReachTrajectory_CurlOrBaseline(tParams);
        tParams.curlGain = -20;
    
        % Store in Array
        AlphaSweep{jj,kk}.Params = tParams;
        AlphaSweep{jj,kk}.BestX = BestX;
%         AlphaSweep{jj,kk}.SOL = SOL_new;
        AlphaSweep{jj,kk}.SOL_LF2 = SOL_LF2_new;
        AlphaSweep{jj,kk}.SOL_EN2 = SOL_EN2_new;
        [AlphaSweep{jj,kk}.costs,AlphaSweep{jj,kk}.costRatio] = calcCostRatios(tParams);
        AlphaSweep{jj,kk}.errorMetrics = calcErrorMetrics(tParams);
%         AlphaSweep{jj,kk}.objfxn = qualityOfSingleTrajectory_noCI(BestX,SOL,SOL_Chan,FitIndices,tParams);
        AlphaSweep{jj,kk}.objfxn = qualityOfLF2EN2Trajectory_noCI(BestX,SOL_LF2,SOL_LF2_Chan,SOL_EN2,SOL_EN2_Chan,FitIndices,tParams);

        % Check if fit is good enough
%         if qualityOfSingleTrajectory_noCI(AlphaSweep{jj,kk}.BestX,SOL,SOL_Chan,FitIndices,tParams) < max_obj_fxn
        if qualityOfLF2EN2Trajectory_noCI(AlphaSweep{jj,kk}.BestX,SOL_LF2,SOL_LF2_Chan,SOL_EN2,SOL_EN2_Chan,FitIndices,tParams) < max_obj_fxn
            kk = kk+1;
        end
    end
end


% Compare Trajectories
figure
for kk = 1:n_sols
    for jj = 1:length(alphas)
%         if AlphaSweep{jj,kk}.objfxn < 0.001
            subplot(1,2,1)
            hold on
            plot(AlphaSweep{jj,kk}.SOL_LF2.x(:,1),AlphaSweep{jj,kk}.SOL_LF2.x(:,2),'color',[0.65,0.65,0.65])
            subplot(1,2,2)
            hold on
            plot(AlphaSweep{jj,kk}.SOL_EN2.x(:,1),AlphaSweep{jj,kk}.SOL_EN2.x(:,2),'color',[0.65,0.65,0.65])
%         end
    end
end
subplot(1,2,1)
plot(SOL_LF2.x(:,1),SOL_LF2.x(:,2),'LineWidth',2,'Color',[0.8,0,0.8])
axis equal
for ii = 1:length(alphas) 
    alpha_str{ii} = num2str(alphas(ii)); 
end
alpha_str{length(alphas) + 1} = 'Original';
legend(alpha_str)
xlim([-0.075, 0.025])
xticks([-0.075:0.025:0.025])
ylim([-0.12 0.12])
yticks([-.1:0.05:.1])

subplot(1,2,2)
plot(SOL_EN2.x(:,1),SOL_EN2.x(:,2),'LineWidth',2,'Color',[0.8,0,0.8])
axis equal
for ii = 1:length(alphas) 
    alpha_str{ii} = num2str(alphas(ii)); 
end
alpha_str{length(alphas) + 1} = 'Original';
legend(alpha_str)
xlim([-0.025, 0.075])
xticks([-0.025:0.025:0.075])
ylim([-0.12 0.12])
yticks([-.1:0.05:.1])

% Compare Cost Ratios
for kk = 1:n_sols
    for jj = 1:length(alphas)
%         if AlphaSweep{jj,kk}.objfxn < 0.001
            cost_ratios(jj,kk) = AlphaSweep{jj,kk}.costRatio;
%         else
%             cost_ratios(jj,kk) = NaN;
%         end
    end
end
% nanmean(log(cost_ratios))
mean(log(cost_ratios))
std(log(cost_ratios))

% Compare Learning
for kk = 1:n_sols
    for jj = 1:length(alphas)
%         if AlphaSweep{jj,kk}.objfxn < 0.001
            estgains(jj,kk) = AlphaSweep{jj,kk}.Params.estGain;
%         else
%             estgains(jj,kk) = NaN;
%         end
    end
end
nanmean(estgains)

