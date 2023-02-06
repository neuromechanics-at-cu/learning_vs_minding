%% Fit Young and Old Trajectories


for ii = 1:length(ages)
    %% Load Data
    if ages{ii} == 'Y'
        load('../../../Data/Trajectory Data/Processed Data/Younger-AveragedTrajectoryData-Individual.mat')
    else
        load('../../../Data/Trajectory Data/Processed Data/Older-AveragedTrajectoryData-Individual.mat')
    end

    % Default: process all subjects
    subj_nums = 1:length(yData.LF2.SubjMeanTrLength);

    for nn = 1:length(subj_nums)
        % Check if trials are valid
        if isnan(yData.LF2.SubjMeanTrLength(subj_nums(nn))) || isnan(yData.EN2.SubjMeanTrLength(subj_nums(nn))) 
            disp(['Subj ',num2str(subj_nums(nn)),' is bad.'])
        else
            %% Sweep for each level of learning
            for jj = 1:length(alphas)
                kk = 1;
                attempt = 0;
                while kk <= n
                    disp(['AgeGp: ',ages{ii},' Subj: ',num2str(subj_nums(nn)),', Attempt: ',num2str(attempt),' Alpha: ',num2str(alphas(jj)),' x',num2str(kk)])
                    %% Define Dynamics
                    clearvars -except ages alphas alpha n ii jj kk yData yDataChan FitIndices ...
                                    subj_nums nn attempt AlphaSweep max_objective_fxn
        
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
                
    
                    % Make Parameter Structure
                    Params.x0 = x0;
                    Params.dt = dt;
                    Params.Asys = Asys;
                    Params.B = Bdis;
                    Params.Fsys = Fsys;
                    Params.Fsys_Chan = Fsys_Chan;
                    
                    Params.curlGain = curlGain;
                    Params.trialLength = 0.6*200;
                    Params.estGain = alphas(jj)*curlGain;
    
                    Params.Q = zeros(size(Asys));
                    Params.R = eye(size(Bdis,2));
                    Params.Phi = eye(size(Asys));
    
                    Params.Mass = Mass;
    
    %                     %initialize BestFitParams Array
    %                     bestFitParams = Params;
    
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
    
    %                     % Update bestFitParams
    %                     bestFitParams = Params;
    
                    FitIndices = logical(FitIndices);
                    numX = sum(FitIndices);
    
                    %%% Number of solutions
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
                        
                        if FitIndices(9)==0 && FitIndices(10)==0
                            [X, fval, flag] = fmincon(@(X) qualityOfModelFit_LF2EN2Only_individual(X,yData,yDataChan,FitIndices,Params,subj_nums(nn)),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
                        else
                            [X, fval, flag] = fmincon(@(X) qualityOfModelFit_AllPhases_individual(X,yData,yDataChan,FitIndices,Params,subj_nums(nn)),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
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
                        if FitIndices(9)==0 && FitIndices(10)==0
                            [X, fval, flag] = fmincon(@(X) qualityOfModelFit_LF2EN2Only_individual(X,yData,yDataChan,FitIndices,Params,subj_nums(nn)),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
                        else               
                            [X, fval, flag] = fmincon(@(X) qualityOfModelFit_AllPhases(X,yData,yDataChan,FitIndices,Params,subj_nums(nn)),Xo,AA,bb,Aeq,beq,lb,ub,nonlcon,options);
                         end
                        Results(m,:) = X';
                        Values(m) = fval;
                        Flags(m) = flag;
                    end
    
                    % Check Extremes
                    Results(2*numsols+1,:) = lb';
                    Results(2*numsols+1,end) = ub(end);
                    if FitIndices(9)==0 && FitIndices(10)==0
                        Values(2*numsols+1) = qualityOfModelFit_LF2EN2Only_individual(lb,yData,yDataChan,FitIndices,Params,subj_nums(nn));
                    else                    
                        Values(2*numsols+1) = qualityOfModelFit_AllPhases_individual(lb,yData,yDataChan,FitIndices,Params,subj_nums(nn));                  
                    end
    
%                     Values

                    % Determine Best Fit from Group
%                     [val, in] = min(Values);
                    [val, in] = nanmin(Values);
                    BestX = Results(in,:);
    
                    if in == 2*numsols + 1
                        pause
                        [val,in] = min(Values(1:end-1));
                        BestX = Results(in,:);
                    end
    
                    tParams = Params;
    
                    %% Generate Model Results for Best Fit
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
                        tParams.estGainLN1 = BestX(count);
                        count = count + 1;
                    end
                    if FitIndices(10)
                        tParams.estGainEF1 = BestX(count);
                        count = count + 1;
                    end
                    if FitIndices(11)
                        tParams.estGainLF2 = BestX(count);
                        count = count + 1;
                    end
                    if FitIndices(12)
                        tParams.estGainEN2 = BestX(count);
                        count = count + 1;
                    end
                    if FitIndices(13)
                        tParams.Q(9,9) = BestX(count);
                        tParams.Q(10,10) = BestX(count);
                        count = count + 1;
                    end
                    if FitIndices(14)
                        tParams.Phi(9,9) = BestX(count);
                        tParams.Phi(10,10) = BestX(count);
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
%                     LN1Params.trialLength = yData.LN1.SubjTrLength(subj_nums(nn));
                    LN1Params.trialLength = floor(yData.LN1.SubjMeanTrLength(subj_nums(nn)));
                    if ~isnan(LN1Params.trialLength)
                        subj_ind_LN1 = find(yData.LN1.SubjMeanTrLength(subj_nums(nn)) == yData.LN1.SubjMeanTrLength(~isnan(yData.LN1.SubjMeanTrLength)),1,'first');
                        forceInit_LN1 = Params.Mass*[yData.LN1.AxSubjNorm(1,subj_ind_LN1); 
                                                    yData.LN1.AySubjNorm(1,subj_ind_LN1)] - ...
                                        LN1Params.curlGain*[0 1;-1 0]*[yData.LN1.VxSubjNorm(1,subj_ind_LN1);
                                                                       yData.LN1.VySubjNorm(1,subj_ind_LN1)];
                        LN1Params.x0 = [yData.LN1.PxSubjNorm(1,subj_ind_LN1);
                                    yData.LN1.PySubjNorm(1,subj_ind_LN1);
                                    yData.LN1.VxSubjNorm(1,subj_ind_LN1);
                                    yData.LN1.VySubjNorm(1,subj_ind_LN1);
                                    forceInit_LN1(1);forceInit_LN1(2);
                                    0;
                                    0.1;
                                    0;
                                    0];
                        SOL_LN1 = generateReachTrajectory_CurlOrBaseline(LN1Params);
                        AccMat = inv(Mass)*([SOL_LN1.x(:,5)';SOL_LN1.x(:,6)'] + LN1Params.curlGain*[0 1;-1 0]*[SOL_LN1.x(:,3)';SOL_LN1.x(:,4)']);
                        SOL_LN1.Ax = AccMat(1,:)';
                        SOL_LN1.Ay = AccMat(2,:)';
                    end
                    % EF1
                    EF1Params = tParams;
                    if FitIndices(10)
                        EF1Params.estGain = tParams.estGainEF1;
                    else
%                         EF1Params.estGain = yData.EF1.estGain(subj_nums(nn));
                        EF1Params.estGain = 0;
                    end
                    EF1Params.curlGain = -20;
                    EF1Params.trialLength = floor(yData.EF1.SubjMeanTrLength(subj_nums(nn)));

                    if ~isnan(EF1Params.trialLength)
                        subj_ind_EF1 = find(yData.EF1.SubjMeanTrLength(subj_nums(nn)) == yData.EF1.SubjMeanTrLength(~isnan(yData.EF1.SubjMeanTrLength)),1,'first');
                        forceInit_EF1 = Params.Mass*[yData.EF1.AxSubjNorm(1,subj_ind_EF1); 
                                                    yData.EF1.AySubjNorm(1,subj_ind_EF1)] - ...
                                        EF1Params.curlGain*[0 1;-1 0]*[yData.EF1.VxSubjNorm(1,subj_ind_EF1);
                                                                       yData.EF1.VySubjNorm(1,subj_ind_EF1)];
                        EF1Params.x0 = [yData.EF1.PxSubjNorm(1,subj_ind_EF1);
                                    yData.EF1.PySubjNorm(1,subj_ind_EF1);
                                    yData.EF1.VxSubjNorm(1,subj_ind_EF1);
                                    yData.EF1.VySubjNorm(1,subj_ind_EF1);
                                    forceInit_EF1(1);forceInit_EF1(2);
                                    0;
                                    0.1;
                                    0;
                                    0];
                        SOL_EF1 = generateReachTrajectory_CurlOrBaseline(EF1Params);
                        AccMat = inv(Mass)*([SOL_EF1.x(:,5)';SOL_EF1.x(:,6)'] + EF1Params.curlGain*[0 1;-1 0]*[SOL_EF1.x(:,3)';SOL_EF1.x(:,4)']);
                        SOL_EF1.Ax = AccMat(1,:)';
                        SOL_EF1.Ay = AccMat(2,:)'; 
                    end
                    % LF2
                    LF2Params = tParams;
                    if FitIndices(11)
                        LF2Params.estGain = tParams.estGainLF2;
                    else
                        LF2Params.estGain = Params.estGain;
                    end
                    LF2Params.curlGain = -20;
                    LF2Params.trialLength = floor(yData.LF2.SubjMeanTrLength(subj_nums(nn)));
                    if ~isnan(LF2Params.trialLength)

                        subj_ind_LF2 = find(yData.LF2.SubjMeanTrLength(subj_nums(nn)) == yData.LF2.SubjMeanTrLength(~isnan(yData.LF2.SubjMeanTrLength)),1,'first');
                        forceInit_LF2 = Params.Mass*[yData.LF2.AxSubjNorm(1,subj_ind_LF2); 
                                                    yData.LF2.AySubjNorm(1,subj_ind_LF2)] - ...
                                        LF2Params.curlGain*[0 1;-1 0]*[yData.LF2.VxSubjNorm(1,subj_ind_LF2);
                                                                       yData.LF2.VySubjNorm(1,subj_ind_LF2)];
                        LF2Params.x0 = [yData.LF2.PxSubjNorm(1,subj_ind_LF2);
                                    yData.LF2.PySubjNorm(1,subj_ind_LF2);
                                    yData.LF2.VxSubjNorm(1,subj_ind_LF2);
                                    yData.LF2.VySubjNorm(1,subj_ind_LF2);
                                    forceInit_LF2(1);forceInit_LF2(2);
                                    0;
                                    0.1;
                                    0;
                                    0]; 
                        SOL_LF2 = generateReachTrajectory_CurlOrBaseline(LF2Params);
                        AccMat = inv(Mass)*([SOL_LF2.x(:,5)';SOL_LF2.x(:,6)'] + LF2Params.curlGain*[0 1;-1 0]*[SOL_LF2.x(:,3)';SOL_LF2.x(:,4)']);
                        SOL_LF2.Ax = AccMat(1,:)';
                        SOL_LF2.Ay = AccMat(2,:)';
                    end
                    % LF2 CHANNEL
                    LF2ChanParams = tParams;
                    if FitIndices(11)
                        LF2ChanParams.estGain = tParams.estGainLF2;
                    else
                        LF2ChanParams.estGain = -20*alpha(mm);
                    end
                    if isnan(yDataChan.LF2.SubjMeanTrLength(subj_nums(nn)))
                        subj_ind_LF2_Chan = NaN;
                    else
                        subj_ind_LF2_Chan = sum(~isnan(yDataChan.LF2.SubjMeanTrLength(1:subj_nums(nn))));
                    end

                    if ~isnan(yDataChan.LF2.SubjMeanTrLength(subj_nums(nn)))
                        LF2ChanParams.curlGain = -20;
                        LF2ChanParams.trialLength = find(~isnan(yDataChan.LF2.PxSubjNorm(:,subj_ind_LF2_Chan)),1,'last');
                        LF2ChanParams.x0 = [yDataChan.LF2.PxSubjNorm(1,subj_ind_LF2_Chan);
                                            yDataChan.LF2.PySubjNorm(1,subj_ind_LF2_Chan);
                                            yDataChan.LF2.VxSubjNorm(1,subj_ind_LF2_Chan);
                                            yDataChan.LF2.VySubjNorm(1,subj_ind_LF2_Chan);
                                            yDataChan.LF2.FxSubjNorm(1,subj_ind_LF2_Chan);
                                            yDataChan.LF2.FySubjNorm(1,subj_ind_LF2_Chan);
                                            0;
                                            0.1;
                                            0;
                                            0];
                        SOL_LF2Chan = generateReachTrajectory_Channel(LF2ChanParams);
                        ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_LF2Chan.x(:,1)';SOL_LF2Chan.x(:,2)';SOL_LF2Chan.x(:,3)';SOL_LF2Chan.x(:,4)'];
                        SOL_LF2Chan.Fx = ForceMat(1,:)';
                        SOL_LF2Chan.Fy = ForceMat(2,:)';
                        SOL_LF2Chan.AdaptInd= fitgain(SOL_LF2Chan.Fx,SOL_LF2Chan.x(:,4)); 
                    end

                    % EN2
                    EN2Params = tParams;
                    if FitIndices(12)
                        EN2Params.estGain = tParams.estGainEN2;
                    elseif FitIndices(11)
                        EN2Params.estGain = tParams.estGainLF2;
                    else
                        EN2Params.estGain = yData.EN2.estGain(subj_nums(nn));
                    end
                    EN2Params.curlGain = 0;
                    EN2Params.trialLength = floor(yData.EN2.SubjMeanTrLength(subj_nums(nn)));
                    if isnan(yData.EN2.SubjMeanTrLength(subj_nums(nn)))
                        subj_ind_EN2 = NaN;
                    else
                        subj_ind_EN2 = sum(~isnan(yData.EN2.SubjMeanTrLength(1:subj_nums(nn))));
                    end
                    if ~isnan(EN2Params.trialLength)
%                         subj_ind_EN2 = find(yData.EN2.SubjMeanTrLength(subj_nums(nn)) == yData.EN2.SubjMeanTrLength(~isnan(yData.EN2.SubjMeanTrLength)),1,'first');    
                        forceInit_EN2 = Params.Mass*[yData.EN2.AxSubjNorm(1,subj_ind_EN2); 
                                                    yData.EN2.AySubjNorm(1,subj_ind_EN2)] - ...
                                        EN2Params.curlGain*[0 1;-1 0]*[yData.EN2.VxSubjNorm(1,subj_ind_EN2);
                                                                    yData.EN2.VySubjNorm(1,subj_ind_EN2)];
                        EN2Params.x0 = [yData.EN2.PxSubjNorm(1,subj_ind_EN2);
                                    yData.EN2.PySubjNorm(1,subj_ind_EN2);
                                    yData.EN2.VxSubjNorm(1,subj_ind_EN2);
                                    yData.EN2.VySubjNorm(1,subj_ind_EN2);
                                    forceInit_EN2(1);forceInit_EN2(2);
                                    0;
                                    0.1;
                                    0;
                                    0];
                        SOL_EN2 = generateReachTrajectory_CurlOrBaseline(EN2Params);
                        AccMat = inv(Mass)*([SOL_EN2.x(:,5)';SOL_EN2.x(:,6)'] + EN2Params.curlGain*[0 1;-1 0]*[SOL_EN2.x(:,3)';SOL_EN2.x(:,4)']);
                        SOL_EN2.Ax = AccMat(1,:)';
                        SOL_EN2.Ay = AccMat(2,:)';
                    end
                    % EN2 CHANNEL
                    EN2ChanParams = tParams;
                    if FitIndices(12)
                        EN2ChanParams.estGain = tParams.estGainEN2;
                    elseif FitIndices(11)
                        EN2ChanParams.estGain = tParams.estGainLF2;
                    else
                        EN2ChanParams.estGain = yData.EN2.estGainMean;
                    end
                    if isnan(yDataChan.EN2.SubjMeanTrLength(subj_nums(nn)))
                        subj_ind_EN2_Chan = NaN;
                    else
                        subj_ind_EN2_Chan = sum(~isnan(yDataChan.EN2.SubjMeanTrLength(1:subj_nums(nn))));
                    end
                    if ~isnan(yDataChan.EN2.SubjMeanTrLength(subj_nums(nn)))
                        EN2ChanParams.curlGain = 0;
                        EN2ChanParams.trialLength = find(~isnan(yDataChan.EN2.PxSubjNorm(:,subj_ind_EN2_Chan)),1,'last');
                        EN2ChanParams.x0 = [yDataChan.EN2.PxSubjNorm(1,subj_ind_EN2_Chan);
                                        yDataChan.EN2.PySubjNorm(1,subj_ind_EN2_Chan);
                                        yDataChan.EN2.VxSubjNorm(1,subj_ind_EN2_Chan);
                                        yDataChan.EN2.VySubjNorm(1,subj_ind_EN2_Chan);
                                        yDataChan.EN2.FxSubjNorm(1,subj_ind_EN2_Chan);
                                        yDataChan.EN2.FySubjNorm(1,subj_ind_EN2_Chan);
                                        0;
                                        0.1;
                                        0;
                                        0];
                        SOL_EN2Chan = generateReachTrajectory_Channel(EN2ChanParams);
                        ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_EN2Chan.x(:,1)';SOL_EN2Chan.x(:,2)';SOL_EN2Chan.x(:,3)';SOL_EN2Chan.x(:,4)'];
                        SOL_EN2Chan.Fx = ForceMat(1,:)';
                        SOL_EN2Chan.Fy = ForceMat(2,:)';
                        SOL_EN2Chan.AdaptInd= fitgain(SOL_EN2Chan.Fx,SOL_EN2Chan.x(:,4));     
                    end
                    % 
                    AlphaSweep{jj,kk}.Params = tParams;
                    AlphaSweep{jj,kk}.BestX = BestX;
                    AlphaSweep{jj,kk}.ObjFxn = min(Values);
                    if ~isnan(LN1Params.trialLength)
                        AlphaSweep{jj,kk}.LN1.SOL = SOL_LN1;
                        AlphaSweep{jj,kk}.LN1.errorMetrics = calcErrorMetrics(LN1Params);
                    end
                    if ~isnan(EF1Params.trialLength)
                        AlphaSweep{jj,kk}.EF1.SOL = SOL_EF1;
                        AlphaSweep{jj,kk}.EF1.errorMetrics = calcErrorMetrics(EF1Params);
                    end
                    if ~isnan(LF2Params.trialLength)
                        AlphaSweep{jj,kk}.LF2.SOL = SOL_LF2;
                        AlphaSweep{jj,kk}.LF2.errorMetrics = calcErrorMetrics(LF2Params);
                        AlphaSweep{jj,kk}.costRatio = calcCostRatios(LF2Params);
                    end
                    if ~isnan(EN2Params.trialLength)
                        AlphaSweep{jj,kk}.EN2.SOL = SOL_EN2;
                        AlphaSweep{jj,kk}.EN2.errorMetrics = calcErrorMetrics(EN2Params);
                    end
                    if ~isnan(yDataChan.LF2.SubjMeanTrLength(subj_nums(nn)))
                        AlphaSweep{jj,kk}.LF2Chan.SOL = SOL_LF2Chan;
                        AlphaSweep{jj,kk}.LF2Chan.errorMetrics = calcErrorMetrics(LF2ChanParams);
                    end
                    if ~isnan(yDataChan.EN2.SubjMeanTrLength(subj_nums(nn)))
                        AlphaSweep{jj,kk}.EN2Chan.SOL = SOL_EN2Chan;
                        AlphaSweep{jj,kk}.EN2Chan.errorMetrics = calcErrorMetrics(EN2ChanParams);
                    end
                    % Check if fit is good enough
                    if min(Values) < max_objective_fxn
                        AlphaSweep{jj,kk}.pass_or_fail = 'pass';
                        kk = kk + 1;
                    else
                        attempt = attempt + 1;
                    end
                    if attempt > 10
                        AlphaSweep{jj,kk}.pass_or_fail = 'fail';
                        kk = kk + 1;
                    end
                end
            end
            save([ages{ii},'-subj',num2str(subj_nums(nn)),'-Sweep',num2str(alphas(1)),'to',num2str(alphas(end)),'-',num2str(now),'.mat'])                    
    
        end
    end
end

