%% Fit Young and Old Trajectories


for ii = 1:length(ages)
    %% Load Data
    if ages{ii} == 'Y'
        load('../../../Data/Trajectory Data/Processed Data/Younger-AveragedTrajectoryData.mat')
    else
        load('../../../Data/Trajectory Data/Processed Data/Older-AveragedTrajectoryData.mat')
    end

    %% Initialize cost weight matrices for each subject group
    q11 = NaN*zeros(length(alphas),n);
    q22 = NaN*zeros(length(alphas),n);
    q3344 = NaN*zeros(length(alphas),n);
    q5566 = NaN*zeros(length(alphas),n);
    q991010 = NaN*zeros(length(alphas),n);
    r1122 = NaN*zeros(length(alphas),n);
    phi1122 = NaN*zeros(length(alphas),n);
    phi3344 = NaN*zeros(length(alphas),n);
    phi5566 = NaN*zeros(length(alphas),n);
    phi991010 = NaN*zeros(length(alphas),n);
    
    %% Sweep for each level of learning
    for jj = 1:length(alphas)

        alpha = alphas(jj)*ones(n,1);
        kk = 1;
        while kk <= n
            alpha(kk)
            kk
            %% Define Dynamics
            clearvars -except ages alphas alpha n ii jj kk q11 q22 ...
                            q3344 q5566 q991010 r1122 phi1122 phi3344 ...
                            phi5566 phi991010 yData yDataChan FitIndices ...
                            ci_scale ci_to_test n_ci

            % Dynamics
            x0 = [0;-0.1;0;0;0;0;0;0.1;0;0];
            curlGain = -20;%Ns/m (curl gain)
            % From experimental Data
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
            
            for mm = 1:length(alpha)

                % Make Parameter Structure
                Params.x0 = x0;
                Params.dt = dt;
                Params.Asys = Asys;
                Params.B = Bdis;
                Params.Fsys = Fsys;
                Params.Fsys_Chan = Fsys_Chan;
                
                Params.curlGain = curlGain;
                Params.trialLength = 0.6*200;
                Params.estGain = alpha(mm)*curlGain;

                Params.Q = zeros(size(Asys));
                Params.R = eye(size(Bdis,2));
                Params.Phi = eye(size(Asys));

                Params.Mass = Mass;

                %initialize BestFitParams Array
                bestFitParams = Params;

                % Initial FMINCON Parameters
                %       [R(1,1), Q(5,5), Q(1,1), Q(2,2), Q(3,3), Phi(1,1), Phi(3,3), Phi(5,5),
                BaseX = [1e2;   0;      1e-10;      0;   0;      1e2;      0;      0;      ...
                        %estGain-LN1,estGain-EF1,estGain-LF2,estGain-EN2;  rate of Force; terminal Rate of Force]
                        0;           0;           Params.estGain;   Params.estGain;         0;              0 ]; 

                %Assign Base Values to Initial Params
                Params.Q(1,1) = BaseX(3); Params.Q(7,7) = BaseX(3);
                Params.Q(1,7) = -1*BaseX(3); Params.Q(7,1) = -1*BaseX(3);
                Params.Q(2,2) = BaseX(4); Params.Q(8,8) = BaseX(4);
                Params.Q(2,8) = -1*BaseX(4); Params.Q(8,2) = -1*BaseX(4);
                Params.Q(3,3) = BaseX(5);
                Params.Q(5,5) = BaseX(2); Params.Q(6,6) = BaseX(2);
                Params.Q(9,9) = BaseX(13); Params.Q(10,10) = BaseX(13);

                Params.R(1,1) = BaseX(1); Params.R(2,2) = BaseX(1);

                Params.Phi(1,1) = BaseX(6); Params.Phi(7,7) = BaseX(6);
                Params.Phi(7,1) = -1*BaseX(6); Params.Phi(1,7) = -1*BaseX(6);
                Params.Phi(2,2) = BaseX(6); Params.Phi(8,8) = BaseX(6);
                Params.Phi(2,8) = -1*BaseX(6); Params.Phi(8,2) = -1*BaseX(6);
                Params.Phi(3,3) = BaseX(7); Params.Phi(4,4) = BaseX(7);
                Params.Phi(5,5) = BaseX(8); Params.Phi(6,6) = BaseX(8);
                Params.Phi(9,9) = BaseX(14); Params.Phi(10,10) = BaseX(14);

                % Update bestFitParams
                bestFitParams = Params;

                FitIndices = logical(FitIndices);
                numX = sum(FitIndices);

                %%% Number of solutions
%                 numsols = 40;
%                 numsols = 20;
                numsols = 10;


                AA = [];
                bb = [];
                Aeq = [];
                beq = [];
                lb = zeros(length(FitIndices),1);
                lb(1) = 1;
                lb(9) = 0;
                lb(10:12) = -20;
                lb = lb(FitIndices);
                ub = 1e20*ones(length(FitIndices),1);
                ub(9:12) = 0;
                ub = ub(FitIndices);
                nonlcon = [];
                options = optimset('MaxFunEval',1E5,'MaxIter',1E5,'TolFun',1E-14,...
                        'TolX',1E-14,'TolCon',1E-14,'LargeScale','on','Algorithm','sqp');
%                 options = optimset('MaxFunEval',1E5,'MaxIter',1E5,'TolFun',1E-16,...
%                         'TolX',1E-16,'TolCon',1E-16,'LargeScale','on','Algorithm','sqp');
                    
                % Initialize Variables
                Values = zeros(numsols,1);
                Results = zeros(numsols,numX);
                Flags = zeros(numsols,1);

                for m = 1:numsols 
                    disp(['minimization loop: ',num2str(m)])

                    Xo = BaseX;
                    % Randomize Starting Point
                    Xo(1:8) = Xo(1:8) + 10*rand(8,1).*Xo(1:8);
                    Xo(9:12) = Xo(9:12) + 2*(rand(4,1)-0.5).*Xo(9:12);
                    Xo(13:14) = Xo(13:14) + 10*rand(2,1).*Xo(13:14);

                    Xo = Xo(FitIndices);
                    
                    if FitIndices(9)==0 && FitIndices(10)==0 && FitIndices(11)==0
                        [X, fval, flag] = fmincon(@(X) qualityOfModelFit_LF2EN2Only(X,yData,yDataChan,FitIndices,Params,ci_scale),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
                    else
                        % When fitting all for phases, we explored three
                        % options for match searching:
                        % 1: Using only error metrics
                        [X, fval, flag] = fmincon(@(X) qualityOfModelFit_AllPhases(X,yData,yDataChan,FitIndices,Params,ci_scale),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
                        % 2: Using position vs time matching
%                         [X, fval, flag] = fmincon(@(X) qualityOfModelFit_AllPhases_Position(X,yData,yDataChan,FitIndices,Params),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
                        % 3: Using position and velocity vs time matching
%                         [X, fval, flag] = fmincon(@(X) qualityOfModelFit_AllPhases_PosVel(X,yData,yDataChan,FitIndices,Params),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);

                    end
                    Results(m,:) = X';
                    Values(m) = fval;
                    Flags(m) = flag;
                end

                % Refined Search
                options = optimset('MaxFunEval',1E4,'MaxIter',1E4,'TolFun',1E-18,...
                    'TolX',1E-16,'TolCon',1E-16,'LargeScale','on','Algorithm','sqp');
                
                for m = numsols+1:2*numsols
                    disp(['refining minimization loop: ',num2str(m)])

                    Xo = Results(m-numsols,:)';
                    if FitIndices(9)==0 && FitIndices(10)==0 && FitIndices(11)==0
                        [X, fval, flag] = fmincon(@(X) qualityOfModelFit_LF2EN2Only(X,yData,yDataChan,FitIndices,Params,ci_scale),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
                    else
                        % When fitting all for phases, we explored three
                        % options for match searching:
                        % 1: Using only error metrics                        
                        [X, fval, flag] = fmincon(@(X) qualityOfModelFit_AllPhases(X,yData,yDataChan,FitIndices,Params,ci_scale),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
                        % 2: Using position vs time matching
%                         [X, fval, flag] = fmincon(@(X) qualityOfModelFit_AllPhases_Position(X,yData,yDataChan,FitIndices,Params),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
                        % 3: Using position and velocity vs time matching
%                         [X, fval, flag] = fmincon(@(X) qualityOfModelFit_AllPhases_PosVel(X,yData,yDataChan,FitIndices,Params),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
                    end
                    Results(m,:) = X';
                    Values(m) = fval;
                    Flags(m) = flag;
                end

                Results(2*numsols+1,:) = lb';
                Results(2*numsols+1,end) = ub(end);
                if FitIndices(9)==0 && FitIndices(10)==0 && FitIndices(11)==0
                    Values(2*numsols+1) = qualityOfModelFit_LF2EN2Only(lb,yData,yDataChan,FitIndices,Params,ci_scale);
                else
                    % When fitting all for phases, we explored three
                    % options for match searching:
                    % 1: Using only error metrics                        
                    Values(2*numsols+1) = qualityOfModelFit_AllPhases(lb,yData,yDataChan,FitIndices,Params,ci_scale);
                    % 2: Using position vs time matching
%                     Values(2*numsols+1) = qualityOfModelFit_AllPhases_Position(lb,yData,yDataChan,FitIndices,Params);                    
                    % 3: Using position and velocity vs time matching
%                     Values(2*numsols+1) = qualityOfModelFit_AllPhases_PosVel(lb,yData,yDataChan,FitIndices,Params);                    
                end
                [val, in] = min(Values);
                BestX = Results(in,:);

                if in == 2*numsols + 1
                    pause
                    [val,in] = min(Values(1:end-1));
                    BestX = Results(in,:);
                end

                tParams = Params;
                for m = 1:size(Results,1)
                    count = 1;
                    if FitIndices(1)
                        tParams.R(1,1) = Results(m,count);
                        tParams.R(2,2) = Results(m,count);
                        count = count + 1;
                    end
                    if FitIndices(2)
                        tParams.Q(5,5) = Results(m,count);
                        tParams.Q(6,6) = Results(m,count);
                        count = count + 1;
                    end
                    if FitIndices(3)
                        tParams.Q(1,1) = Results(m,count);
                        tParams.Q(7,7) = Results(m,count);
                        tParams.Q(1,7) = -1*Results(m,count);
                        tParams.Q(7,1) = -1*Results(m,count);
                        count = count + 1;
                    end
                    if FitIndices(4) 
                        tParams.Q(2,2) = Results(m,count);
                        tParams.Q(8,8) = Results(m,count);
                        tParams.Q(2,8) = -1*Results(m,count);
                        tParams.Q(8,2) = -1*Results(m,count);
                        count = count + 1;
                    end
                    if FitIndices(5)
                       tParams.Q(3,3) = Results(m,count); 
                       count = count + 1;
                    end
                    if FitIndices(6)
                        tParams.Phi(1,1) = Results(m,count);
                        tParams.Phi(7,7) = Results(m,count);
                        tParams.Phi(1,7) = -1*Results(m,count);
                        tParams.Phi(7,1) = -1*Results(m,count);
                        tParams.Phi(2,2) = Results(m,count);
                        tParams.Phi(8,8) = Results(m,count);
                        tParams.Phi(2,8) = -1*Results(m,count);
                        tParams.Phi(8,2) = -1*Results(m,count);
                        count = count + 1;
                    end
                    if FitIndices(7)
                        tParams.Phi(3,3) = Results(m,count);
                        tParams.Phi(4,4) = Results(m,count);
                        count = count + 1;
                    end
                    if FitIndices(8)
                        tParams.Phi(5,5) = Results(m,count);
                        tParams.Phi(6,6) = Results(m,count);
                        count = count + 1;
                    end
                    if FitIndices(9)
                        tParams.estGainLN1 = Results(m,count);
                        count = count + 1;
                    end
                    if FitIndices(10)
                        tParams.estGainEF1 = Results(m,count);
                        count = count + 1;
                    end
                    if FitIndices(11)
                        tParams.estGainLF2 = Results(m,count);
                        count = count + 1;
                    end
                    if FitIndices(12)
                        tParams.estGainEN2 = Results(m,count);
                        count = count + 1;
                    end
                    if FitIndices(13)
                        tParams.Q(9,9) = Results(m,count);
                        tParams.Q(10,10) = Results(m,count);
                        count = count + 1;
                    end
                    if FitIndices(14)
                        tParams.Phi(9,9) = Results(m,count);
                        tParams.Phi(10,10) = Results(m,count);
                    end
    
                    % Generate Best Fit Trajectories
                    % LN1
                    LN1Params = tParams;
                    if FitIndices(9)
                        LN1Params.estGain = tParams.estGainLN1;
                    else
                        LN1Params.estGain = 0;
                    end
                    LN1Params.curlGain = 0;
                    LN1Params.trialLength = length(yData.LN1.PxSubjMean);
                    forceInit_LN1 = Params.Mass*[yData.LN1.AxSubjMean(1); 
                                                yData.LN1.AySubjMean(1)] - ...
                                    LN1Params.curlGain*[0 1;-1 0]*[yData.LN1.VxSubjMean(1);
                                                                   yData.LN1.VySubjMean(1)];
                    LN1Params.x0 = [yData.LN1.PxSubjMean(1);
                                yData.LN1.PySubjMean(1);
                                yData.LN1.VxSubjMean(1);
                                yData.LN1.VySubjMean(1);
                                forceInit_LN1(1);forceInit_LN1(2);
                                0;
                                0.1;
                                0;
                                0];
                    SOL_LN1 = generateReachTrajectory_CurlOrBaseline(LN1Params);
                    AccMat = inv(Mass)*([SOL_LN1.x(:,5)';SOL_LN1.x(:,6)'] + LN1Params.curlGain*[0 1;-1 0]*[SOL_LN1.x(:,3)';SOL_LN1.x(:,4)']);
                    SOL_LN1.Ax = AccMat(1,:)';
                    SOL_LN1.Ay = AccMat(2,:)';
                    Young.LN1.AllSol.Params(m) = LN1Params;
                    Young.LN1.AllSol.Traj(m) = SOL_LN1;
                    % EF1
                    EF1Params = tParams;
                    if FitIndices(10)
                        EF1Params.estGain = tParams.estGainEF1;
                    else
                        EF1Params.estGain = yData.EF1.estGainMean;
                    end
                    EF1Params.curlGain = -20;
                    EF1Params.trialLength = length(yData.EF1.PxSubjMean);
                    forceInit_EF1 = Params.Mass*[yData.EF1.AxSubjMean(1); 
                                                yData.EF1.AySubjMean(1)] - ...
                                    EF1Params.curlGain*[0 1;-1 0]*[yData.EF1.VxSubjMean(1);
                                                                   yData.EF1.VySubjMean(1)];
                    EF1Params.x0 = [yData.EF1.PxSubjMean(1);
                                yData.EF1.PySubjMean(1);
                                yData.EF1.VxSubjMean(1);
                                yData.EF1.VySubjMean(1);
                                forceInit_EF1(1);forceInit_EF1(2);
                                0;
                                0.1;
                                0;
                                0];
                    SOL_EF1 = generateReachTrajectory_CurlOrBaseline(EF1Params);
                    AccMat = inv(Mass)*([SOL_EF1.x(:,5)';SOL_EF1.x(:,6)'] + EF1Params.curlGain*[0 1;-1 0]*[SOL_EF1.x(:,3)';SOL_EF1.x(:,4)']);
                    SOL_EF1.Ax = AccMat(1,:)';
                    SOL_EF1.Ay = AccMat(2,:)';
                    Young.EF1.AllSol.Params(m) = EF1Params;
                    Young.EF1.AllSol.Traj(m) = SOL_EF1;
                    
                    % LF2
                    LF2Params = tParams;
                    if FitIndices(11)
                        LF2Params.estGain = tParams.estGainLF2;
                    else
                        LF2Params.estGain = Params.estGain;
                    end
                    LF2Params.curlGain = -20;
                    LF2Params.trialLength = length(yData.LF2.PxSubjMean);
                    forceInit_LF2 = Params.Mass*[yData.LF2.AxSubjMean(1); 
                                                yData.LF2.AySubjMean(1)] - ...
                                    LF2Params.curlGain*[0 1;-1 0]*[yData.LF2.VxSubjMean(1);
                                                                   yData.LF2.VySubjMean(1)];
                    LF2Params.x0 = [yData.LF2.PxSubjMean(1);
                                yData.LF2.PySubjMean(1);
                                yData.LF2.VxSubjMean(1);
                                yData.LF2.VySubjMean(1);
                                forceInit_LF2(1);forceInit_LF2(2);
                                0;
                                0.1;
                                0;
                                0]; 
                    SOL_LF2 = generateReachTrajectory_CurlOrBaseline(LF2Params);
                    AccMat = inv(Mass)*([SOL_LF2.x(:,5)';SOL_LF2.x(:,6)'] + LF2Params.curlGain*[0 1;-1 0]*[SOL_LF2.x(:,3)';SOL_LF2.x(:,4)']);
                    SOL_LF2.Ax = AccMat(1,:)';
                    SOL_LF2.Ay = AccMat(2,:)';
                    Young.LF2.AllSol.Params(m) = LF2Params;
                    Young.LF2.AllSol.Traj(m) = SOL_LF2;
    
                    % EN2
                    EN2Params = tParams;
                    if FitIndices(12)
                        EN2Params.estGain = tParams.estGainEN2;
                    else
                        EN2Params.estGain = yData.EN2.estGainMean;
                    end
                    EN2Params.curlGain = 0;
                    EN2Params.trialLength = length(yData.EN2.PxSubjMean);
                    forceInit_EN2 = Params.Mass*[yData.EN2.AxSubjMean(1); 
                                                yData.EN2.AySubjMean(1)] - ...
                                    EN2Params.curlGain*[0 1;-1 0]*[yData.EN2.VxSubjMean(1);
                                                                yData.EN2.VySubjMean(1)];
                    EN2Params.x0 = [yData.EN2.PxSubjMean(1);
                                yData.EN2.PySubjMean(1);
                                yData.EN2.VxSubjMean(1);
                                yData.EN2.VySubjMean(1);
                                forceInit_EN2(1);forceInit_EN2(2);
                                0;
                                0.1;
                                0;
                                0];
                    SOL_EN2 = generateReachTrajectory_CurlOrBaseline(EN2Params);
                    AccMat = inv(Mass)*([SOL_EN2.x(:,5)';SOL_EN2.x(:,6)'] + EN2Params.curlGain*[0 1;-1 0]*[SOL_EN2.x(:,3)';SOL_EN2.x(:,4)']);
                    SOL_EN2.Ax = AccMat(1,:)';
                    SOL_EN2.Ay = AccMat(2,:)';
                    Young.EN2.AllSol.Params(m) = EN2Params;
                    Young.EN2.AllSol.Traj(m) = SOL_EN2;
                end
                
                count = 1;
                if FitIndices(1)
                    bestFitParams.R(1,1) = BestX(count);
                    bestFitParams.R(2,2) = BestX(count);
                    count = count + 1;
                end
                if FitIndices(2)
                    bestFitParams.Q(5,5) = BestX(count);
                    bestFitParams.Q(6,6) = BestX(count);
                    count = count + 1;
                end
                if FitIndices(3)
                    bestFitParams.Q(1,1) = BestX(count);
                    bestFitParams.Q(7,7) = BestX(count);
                    bestFitParams.Q(1,7) = -1*BestX(count);
                    bestFitParams.Q(7,1) = -1*BestX(count);
                    count = count + 1;
                end
                if FitIndices(4) 
                    bestFitParams.Q(2,2) = BestX(count);
                    bestFitParams.Q(8,8) = BestX(count);
                    bestFitParams.Q(2,8) = -1*BestX(count);
                    bestFitParams.Q(8,2) = -1*BestX(count);
                    count = count + 1;
                end
                if FitIndices(5)
                   bestFitParams.Q(3,3) = BestX(count); 
                   count = count + 1;
                end
                if FitIndices(6)
                    bestFitParams.Phi(1,1) = BestX(count);
                    bestFitParams.Phi(7,7) = BestX(count);
                    bestFitParams.Phi(1,7) = -1*BestX(count);
                    bestFitParams.Phi(7,1) = -1*BestX(count);
                    bestFitParams.Phi(2,2) = BestX(count);
                    bestFitParams.Phi(8,8) = BestX(count);
                    bestFitParams.Phi(2,8) = -1*BestX(count);
                    bestFitParams.Phi(8,2) = -1*BestX(count);
                    count = count + 1;
                end
                if FitIndices(7)
                    bestFitParams.Phi(3,3) = BestX(count);
                    bestFitParams.Phi(4,4) = BestX(count);
                    count = count + 1;
                end
                if FitIndices(8)
                    bestFitParams.Phi(5,5) = BestX(count);
                    bestFitParams.Phi(6,6) = BestX(count);
                    count = count + 1;
                end
                if FitIndices(9)
                    bestFitParams.estGainLN1 = BestX(count);
                    count = count + 1;
                end
                if FitIndices(10)
                    bestFitParams.estGainEF1 = BestX(count);
                    count = count + 1;
                end
                if FitIndices(11)
                    bestFitParams.estGainLF2 = BestX(count);
                    count = count + 1;
                end
                if FitIndices(12)
                    bestFitParams.estGainEN2 = BestX(count);
                    count = count + 1;
                end
                if FitIndices(13)
                    bestFitParams.Q(9,9) = BestX(count);
                    bestFitParams.Q(10,10) = BestX(count);
                    count = count + 1;
                end
                if FitIndices(14)
                    bestFitParams.Phi(9,9) = BestX(count);
                    bestFitParams.Phi(10,10) = BestX(count);
                end


                % Generate Best Fit Trajectories
                LN1Params = bestFitParams;
                if FitIndices(9)
                    LN1Params.estGain = bestFitParams.estGainLN1;
                else
                    LN1Params.estGain = 0;
                end
                LN1Params.curlGain = 0;
                LN1Params.trialLength = length(yData.LN1.PxSubjMean);
                forceInit_LN1 = Params.Mass*[yData.LN1.AxSubjMean(1); 
                                            yData.LN1.AySubjMean(1)] - ...
                                LN1Params.curlGain*[0 1;-1 0]*[yData.LN1.VxSubjMean(1);
                                                            yData.LN1.VySubjMean(1)];
                LN1Params.x0 = [yData.LN1.PxSubjMean(1);
                            yData.LN1.PySubjMean(1);
                            yData.LN1.VxSubjMean(1);
                            yData.LN1.VySubjMean(1);
                            forceInit_LN1(1);forceInit_LN1(2);
                            0;
                            0.1;
                            0;
                            0]; 
                SOL_LN1 = generateReachTrajectory_CurlOrBaseline(LN1Params);
                AccMat = inv(Mass)*([SOL_LN1.x(:,5)';SOL_LN1.x(:,6)'] + LN1Params.curlGain*[0 1;-1 0]*[SOL_LN1.x(:,3)';SOL_LN1.x(:,4)']);
                SOL_LN1.Ax = AccMat(1,:)';
                SOL_LN1.Ay = AccMat(2,:)';
                Young.LN1.Params = LN1Params;
                Young.LN1.Traj = SOL_LN1;
                % LN1 Channel
                LN1ChanParams = bestFitParams;
                if FitIndices(9)
                    LN1ChanParams.estGain = bestFitParams.estGainLN1;
                else
                    LN1ChanParams.estGain = 0;
                end
                LN1ChanParams.curlGain = 0;
                LN1ChanParams.trialLength = length(yDataChan.LN1.PxSubjMean);
                LN1ChanParams.curlGain = 0;
                LN1ChanParams.trialLength = length(yDataChan.LN1.PxSubjMean);
                LN1ChanParams.x0 = [yDataChan.LN1.PxSubjMean(1);
                                    yDataChan.LN1.PySubjMean(1);
                                    yDataChan.LN1.VxSubjMean(1);
                                    yDataChan.LN1.VySubjMean(1);
                                    yDataChan.LN1.FxSubjMean(1);
                                    yDataChan.LN1.FySubjMean(1);
                                    0;
                                    0.1;0;0];
                SOL_LN1Chan = generateReachTrajectory_Channel(LN1ChanParams);
                ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_LN1Chan.x(:,1)';SOL_LN1Chan.x(:,2)';SOL_LN1Chan.x(:,3)';SOL_LN1Chan.x(:,4)'];
                SOL_LN1Chan.Fx = ForceMat(1,:)';
                SOL_LN1Chan.Fy = ForceMat(2,:)';
                SOL_LN1Chan.AdaptInd= fitgain(SOL_LN1Chan.Fx,SOL_LN1Chan.x(:,4)); 
                Young.LN1Chan.Params = LN1ChanParams;
                Young.LN1Chan.Traj = SOL_LN1Chan;

                %EF1
                EF1Params = bestFitParams;
                if FitIndices(10)
                    EF1Params.estGain = bestFitParams.estGainEF1;
                else
                    EF1Params.estGain = yData.EF1.estGainMean;
                end
                EF1Params.curlGain = -20;
                EF1Params.trialLength = length(yData.EF1.PxSubjMean);
                forceInit_EF1 = Params.Mass*[yData.EF1.AxSubjMean(1); 
                                            yData.EF1.AySubjMean(1)] - ...
                                EF1Params.curlGain*[0 1;-1 0]*[yData.EF1.VxSubjMean(1);
                                                            yData.EF1.VySubjMean(1)];
                EF1Params.x0 = [yData.EF1.PxSubjMean(1);
                            yData.EF1.PySubjMean(1);
                            yData.EF1.VxSubjMean(1);
                            yData.EF1.VySubjMean(1);
                            forceInit_EF1(1);forceInit_EF1(2);
                            0;
                            0.1;0;0];
                SOL_EF1 = generateReachTrajectory_CurlOrBaseline(EF1Params);
                AccMat = inv(Mass)*([SOL_EF1.x(:,5)';SOL_EF1.x(:,6)'] + EF1Params.curlGain*[0 1;-1 0]*[SOL_EF1.x(:,3)';SOL_EF1.x(:,4)']);
                SOL_EF1.Ax = AccMat(1,:)';
                SOL_EF1.Ay = AccMat(2,:)';
                Young.EF1.Params = EF1Params;
                Young.EF1.Traj = SOL_EF1;
                % EF1 CHANNEL
                EF1ChanParams = bestFitParams;
                if FitIndices(10)
                    EF1ChanParams.estGain = bestFitParams.estGainEF1;
                else
                    EF1ChanParams.estGain = yData.EF1.estGainMean;
                end
                EF1ChanParams.curlGain = -20;
                EF1ChanParams.trialLength = length(yDataChan.EF1.PxSubjMean);
                EF1ChanParams.x0 = [yDataChan.EF1.PxSubjMean(1);
                                    yDataChan.EF1.PySubjMean(1);
                                    yDataChan.EF1.VxSubjMean(1);
                                    yDataChan.EF1.VySubjMean(1);
                                    yDataChan.EF1.FxSubjMean(1);
                                    yDataChan.EF1.FySubjMean(1);
                                    0;
                                    0.1;0;0];
                SOL_EF1Chan = generateReachTrajectory_Channel(EF1ChanParams);
                ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_EF1Chan.x(:,1)';SOL_EF1Chan.x(:,2)';SOL_EF1Chan.x(:,3)';SOL_EF1Chan.x(:,4)'];
                SOL_EF1Chan.Fx = ForceMat(1,:)';
                SOL_EF1Chan.Fy = ForceMat(2,:)';
                SOL_EF1Chan.AdaptInd= fitgain(SOL_EF1Chan.Fx,SOL_EF1Chan.x(:,4)); 
                Young.EF1Chan.Params = EF1ChanParams;
                Young.EF1Chan.Traj = SOL_EF1Chan;
                
                %LF2
                LF2Params = bestFitParams;
                if FitIndices(11)
                    LF2Params.estGain = bestFitParams.estGainLF2;
                else
                    LF2Params.estGain = Params.estGain;
                end
                LF2Params.curlGain = -20;
                LF2Params.trialLength = length(yData.LF2.PxSubjMean);
                forceInit_LF2 = Params.Mass*[yData.LF2.AxSubjMean(1); 
                                            yData.LF2.AySubjMean(1)] - ...
                                LF2Params.curlGain*[0 1;-1 0]*[yData.LF2.VxSubjMean(1);
                                                            yData.LF2.VySubjMean(1)];
                LF2Params.x0 = [yData.LF2.PxSubjMean(1);
                            yData.LF2.PySubjMean(1);
                            yData.LF2.VxSubjMean(1);
                            yData.LF2.VySubjMean(1);
                            forceInit_LF2(1);forceInit_LF2(2);
                            0;
                            0.1;
                            0;
                            0]; 
                SOL_LF2 = generateReachTrajectory_CurlOrBaseline(LF2Params);
                AccMat = inv(Mass)*([SOL_LF2.x(:,5)';SOL_LF2.x(:,6)'] + LF2Params.curlGain*[0 1;-1 0]*[SOL_LF2.x(:,3)';SOL_LF2.x(:,4)']);
                SOL_LF2.Ax = AccMat(1,:)';
                SOL_LF2.Ay = AccMat(2,:)';
                Young.LF2.Params = LF2Params;
                Young.LF2.Traj = SOL_LF2;
                % LF2 CHANNEL
                LF2ChanParams = bestFitParams;
                if FitIndices(11)
                    LF2ChanParams.estGain = bestFitParams.estGainLF2;
                else
                    LF2ChanParams.estGain = -20*alpha(mm);
                end
                LF2ChanParams.curlGain = -20;
                LF2ChanParams.trialLength = length(yDataChan.LF2.PxSubjMean);
                LF2ChanParams.x0 = [yDataChan.LF2.PxSubjMean(1);
                                    yDataChan.LF2.PySubjMean(1);
                                    yDataChan.LF2.VxSubjMean(1);
                                    yDataChan.LF2.VySubjMean(1);
                                    yDataChan.LF2.FxSubjMean(1);
                                    yDataChan.LF2.FySubjMean(1);
                                    0;
                                    0.1;0;0];
                SOL_LF2Chan = generateReachTrajectory_Channel(LF2ChanParams);
                ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_LF2Chan.x(:,1)';SOL_LF2Chan.x(:,2)';SOL_LF2Chan.x(:,3)';SOL_LF2Chan.x(:,4)'];
                SOL_LF2Chan.Fx = ForceMat(1,:)';
                SOL_LF2Chan.Fy = ForceMat(2,:)';
                SOL_LF2Chan.AdaptInd= fitgain(SOL_LF2Chan.Fx,SOL_LF2Chan.x(:,4)); 
                Young.LF2Chan.Params = LF2ChanParams;
                Young.LF2Chan.Traj = SOL_LF2Chan;
                %EN2
                EN2Params = bestFitParams;
                if FitIndices(12)
                    EN2Params.estGain = bestFitParams.estGainEN2;
                else
                    EN2Params.estGain = yData.EN2.estGainMean;
                end
                EN2Params.curlGain = 0;
                EN2Params.trialLength = length(yData.EN2.PxSubjMean);
                forceInit_EN2 = Params.Mass*[yData.EN2.AxSubjMean(1); 
                                            yData.EN2.AySubjMean(1)] - ...
                                EN2Params.curlGain*[0 1;-1 0]*[yData.EN2.VxSubjMean(1);
                                                            yData.EN2.VySubjMean(1)];
                EN2Params.x0 = [yData.EN2.PxSubjMean(1);
                            yData.EN2.PySubjMean(1);
                            yData.EN2.VxSubjMean(1);
                            yData.EN2.VySubjMean(1);
                            forceInit_EN2(1);forceInit_EN2(2);
                            0;
                            0.1;
                            0;
                            0];
                SOL_EN2 = generateReachTrajectory_CurlOrBaseline(EN2Params);
                AccMat = inv(Mass)*([SOL_EN2.x(:,5)';SOL_EN2.x(:,6)'] + EN2Params.curlGain*[0 1;-1 0]*[SOL_EN2.x(:,3)';SOL_EN2.x(:,4)']);
                SOL_EN2.Ax = AccMat(1,:)';
                SOL_EN2.Ay = AccMat(2,:)';
                Young.EN2.Params = EN2Params;
                Young.EN2.Traj = SOL_EN2;
                % EN2 CHANNEL
                EN2ChanParams = bestFitParams;
                if FitIndices(12)
                    EN2ChanParams.estGain = bestFitParams.estGainEN2;
                else
                    EN2ChanParams.estGain = yData.EN2.estGainMean;
                end
                EN2ChanParams.curlGain = -20;%Not Necessary...
                EN2ChanParams.trialLength = length(yDataChan.EN2.PxSubjMean);
                EN2ChanParams.x0 = [yDataChan.EN2.PxSubjMean(1);
                                    yDataChan.EN2.PySubjMean(1);
                                    yDataChan.EN2.VxSubjMean(1);
                                    yDataChan.EN2.VySubjMean(1);
                                    yDataChan.EN2.FxSubjMean(1);
                                    yDataChan.EN2.FySubjMean(1);
                                    0;
                                    0.1;0;0];
                SOL_EN2Chan = generateReachTrajectory_Channel(EN2ChanParams);
                ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_EN2Chan.x(:,1)';SOL_EN2Chan.x(:,2)';SOL_EN2Chan.x(:,3)';SOL_EN2Chan.x(:,4)'];
                SOL_EN2Chan.Fx = ForceMat(1,:)';
                SOL_EN2Chan.Fy = ForceMat(2,:)';
                SOL_EN2Chan.AdaptInd= fitgain(SOL_EN2Chan.Fx,SOL_EN2Chan.x(:,4)); 
                Young.EN2Chan.Params = EN2ChanParams;
                Young.EN2Chan.Traj = SOL_EN2Chan;

                YoungAlphaSweep(jj,mm) = Young;

                % MAX PERP
                % Of Data
                %LN1
                ind = find(abs(yData.LN1.PxSubjMean) == max(abs(yData.LN1.PxSubjMean)),1,'first');
                maxPerp_LN1 = yData.LN1.PxSubjMean(ind);
                maxPerpSE_LN1 = yData.LN1.PxSubjSE(ind);
                %EF1
                ind = find(abs(yData.EF1.PxSubjMean) == max(abs(yData.EF1.PxSubjMean)),1,'first');
                maxPerp_EF1 = yData.EF1.PxSubjMean(ind);
                maxPerpSE_EF1 = yData.EF1.PxSubjSE(ind);
                %LF2
                ind = find(abs(yData.LF2.PxSubjMean) == max(abs(yData.LF2.PxSubjMean)),1,'first');
                maxPerp_LF2 = yData.LF2.PxSubjMean(ind);
                maxPerpSE_LF2 = yData.LF2.PxSubjSE(ind);
                %EN2
                ind = find(abs(yData.EN2.PxSubjMean) == max(abs(yData.EN2.PxSubjMean)),1,'first');
                maxPerp_EN2 = yData.EN2.PxSubjMean(ind);
                maxPerpSE_EN2 = yData.EN2.PxSubjSE(ind); 
                % Of Alpha Fit
                %LN1
                ind = find(abs(YoungAlphaSweep(jj,mm).LN1.Traj.x(:,1)) == max(abs(YoungAlphaSweep(jj,mm).LN1.Traj.x(:,1))),1,'first');
                maxPerp_mLN1 = YoungAlphaSweep(jj,mm).LN1.Traj.x(ind,1);
                %EF1
                ind = find(abs(YoungAlphaSweep(jj,mm).EF1.Traj.x(:,1)) == max(abs(YoungAlphaSweep(jj,mm).EF1.Traj.x(:,1))),1,'first');
                maxPerp_mEF1 = YoungAlphaSweep(jj,mm).EF1.Traj.x(ind,1);
                %LF2
                ind = find(abs(YoungAlphaSweep(jj,mm).LF2.Traj.x(:,1)) == max(abs(YoungAlphaSweep(jj,mm).LF2.Traj.x(:,1))),1,'first');
                maxPerp_mLF2 = YoungAlphaSweep(jj,mm).LF2.Traj.x(ind,1);
                %EN2
                ind = find(abs(YoungAlphaSweep(jj,mm).EN2.Traj.x(:,1)) == max(abs(YoungAlphaSweep(jj,mm).EN2.Traj.x(:,1))),1,'first');
                maxPerp_mEN2 = YoungAlphaSweep(jj,mm).EN2.Traj.x(ind,1);
    
            end
            %% Stats Analysis
%             if FitIndices(9)==0 && FitIndices(10)==0 && FitIndices(11)==0
%                 phaseNames = {'LF2','EN2'};
%             else
                phaseNames = {'LN1','EF1','LF2','EN2'};                
%             end
            % Max Perp Force
            maxFx = zeros(length(phaseNames),2);
            maxFxSE = zeros(length(phaseNames),2);
            % Adaptation Index
            stiffConst = 2000;
            dampConst = 50;
            adaptInd = zeros(length(phaseNames),2);
            adaptIndSE = zeros(length(phaseNames),2);
            % Max Perp Error
            maxPerp = zeros(length(phaseNames),2);
            maxPerpSE = zeros(length(phaseNames),2);
            
            for mm = 1:length(phaseNames) 
                % MaxFx
                eval(['[~,ind] = max(abs(yDataChan.',phaseNames{mm},'.FxSubjMean));'])
                eval(['maxFx(mm,1) = yDataChan.',phaseNames{mm},'.FxSubjMean(ind);']);
                eval(['maxFxSE(mm,1) = ci_scale*yDataChan.',phaseNames{mm},'.FxSubjSE(ind);']);
                eval(['[~,ind] = max(abs(Young.',phaseNames{mm},'Chan.Traj.Fx));']);
                eval(['maxFx(mm,2) = Young.',phaseNames{mm},'Chan.Traj.Fx(ind);'])
                % Adapt Ind
                eval(['adaptInd(mm,1) = yDataChan.',phaseNames{mm},'.estGainMean;'])
                eval(['adaptIndSE(mm,1) = ci_scale*yDataChan.',phaseNames{mm},'.estGainSE;'])
                eval(['adaptInd(mm,2) = Young.',phaseNames{mm},'Chan.Traj.AdaptInd;'])
                % MaxPx
                eval(['[~,ind] = max(abs(yData.',phaseNames{mm},'.PxSubjMean));']);
                eval(['maxPerp(mm,1) = yData.',phaseNames{mm},'.PxSubjMean(ind);']);
                eval(['maxPerpSE(mm,1) = ci_scale*yData.',phaseNames{mm},'.PxSubjSE(ind);']);
                eval(['[~,ind] = max(abs(Young.',phaseNames{mm},'.Traj.x(:,1)));'])
                eval(['maxPerp(mm,2) = Young.',phaseNames{mm},'.Traj.x(ind,1);'])
            end
            % Make a Binary Matrix of Results
            stats = zeros(1,3,length(phaseNames));
            for mm = 1:length(phaseNames) 
                % MaxFx
                if (maxFx(mm,2) < maxFx(mm,1) + maxFxSE(mm,1)) & (maxFx(mm,2) > maxFx(mm,1) - maxFxSE(mm,1))
                    stats(1,1,mm) = 1;
                end
                % AdaptInd
                if (adaptInd(mm,2) < adaptInd(mm,1) + adaptIndSE(mm,1)) & (adaptInd(mm,2) > adaptInd(mm,1) - adaptIndSE(mm,1))
                    stats(1,2,mm) = 1;
                end    
                % MaxPx
                if (maxPerp(mm,2) < maxPerp(mm,1) + maxPerpSE(mm,1)) & (maxPerp(mm,2) > maxPerp(mm,1) - maxPerpSE(mm,1))
                    stats(1,3,mm) = 1;
                end    
            end
            % AND across Phases and Metrics... see what fits all
            stats_byPhase = squeeze(prod(stats,2));
%             stats_allPhase_allMetrics = prod(stats_byPhase,2);         
            stats_allPhase_allMetrics = prod(stats_byPhase);         
            
            %% Enter Cost Weights into Matrix
            if stats_allPhase_allMetrics
                q11(jj,kk) = Young.LF2.Params.Q(1,1);
                q22(jj,kk) = Young.LF2.Params.Q(2,2);
                q3344(jj,kk) = Young.LF2.Params.Q(3,3);
                q5566(jj,kk) = Young.LF2.Params.Q(5,5);
                q991010(jj,kk) = Young.LF2.Params.Q(9,9);
                r1122(jj,kk) = Young.LF2.Params.R(1,1);
                phi1122(jj,kk) = Young.LF2.Params.Phi(1,1);
                phi3344(jj,kk) = Young.LF2.Params.Phi(3,3);
                phi5566(jj,kk) = Young.LF2.Params.Phi(5,5);
                phi991010(jj,kk) = Young.LF2.Params.Phi(9,9);
                %% Save Result if stats ok
                if FitIndices(9)==0 && FitIndices(10)==0 && FitIndices(11)==0
                    save([ages{ii},'-',num2str(kk),'-alpha',num2str(alpha(kk)),'-',ci_scale,'.mat'])
                else
                    save([ages{ii},'-',num2str(kk),'-',ci_scale,'-fitAllPhases.mat'])                    
                end
                % Since it is a valid solution, increment kk
                kk = kk + 1;
            % In the case of all four phases, LN1 will rarely fit due to it 
            % being a symmetrical point mass so it is still acceptable to 
            % save the fit, if only these phases does not meet specs. We 
            % save the fit if only LF2 and EN2 are acceptable.
%             elseif FitIndices(9)==0 && FitIndices(10)==0 && FitIndices(11)==0 && prod(stats_byPhase(3:end))    
            elseif prod(stats_byPhase(3:end))    
                q11(jj,kk) = Young.LF2.Params.Q(1,1);
                q22(jj,kk) = Young.LF2.Params.Q(2,2);
                q3344(jj,kk) = Young.LF2.Params.Q(3,3);
                q5566(jj,kk) = Young.LF2.Params.Q(5,5);
                q991010(jj,kk) = Young.LF2.Params.Q(9,9);
                r1122(jj,kk) = Young.LF2.Params.R(1,1);
                phi1122(jj,kk) = Young.LF2.Params.Phi(1,1);
                phi3344(jj,kk) = Young.LF2.Params.Phi(3,3);
                phi5566(jj,kk) = Young.LF2.Params.Phi(5,5);
                phi991010(jj,kk) = Young.LF2.Params.Phi(9,9);
                % If you arent fitting over LN1 and EF1, then it could
                % still be a valid solution
                if FitIndices(9)==0 && FitIndices(10)==0 && FitIndices(11)==0
                    save([ages{ii},'-',num2str(kk),'-alpha',num2str(alpha(kk)),'-',ci_scale,'.mat'])
                else
                    % Otherwise, if we're fitting over all phases, we still
                    % want to save the solution, and increment anyway.
                    save([ages{ii},'-',num2str(kk),'-',ci_scale','-fitAllPhases-badLN1.mat'])                    
                end
                % Increment kk
                kk = kk + 1;   
            elseif FitIndices(9)==0 && FitIndices(10)==0 && FitIndices(11)==0 && prod(stats_byPhase(3:end))==0
                save([ages{ii},'-',num2str(kk),'-alpha',num2str(alpha(kk)),'-',ci_scale,'-badfit-',num2str(now),'.mat'])
                % Toggle this on or off, if you'd like to continue on even
                % though the fit is bad.
                kk = kk + 1;
            end
        end
    end
    %% Save Result of All Runs, Per Age, Per Alpha
    save([ages{ii},'-alpha',num2str(alpha),'-CostWeights.mat'],'q11','q22','q3344',...
        'q5566','q991010','r1122','phi1122','phi3344','phi5566','phi991010')
end

