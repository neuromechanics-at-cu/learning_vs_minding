%% Figure 3: Initial Model Fits Across All Phases

if ~exist("ages")
    ages = {'Y','O'};
end

addpath('../../Scripts/Model Fitting')
addpath('../Supporting Functions')


for nn = 1:length(ages);
    clearvars -except ages nn

    %% Load Data
    if ages{nn} == 'Y'
        load('../../Data/Model Fits/Model Accuracy Fits (Fig 4)/BestFitModel-FourPhases-Younger.mat')
        color = [0,128/255,0];
        titletext = 'Younger';
    elseif ages{nn} == 'O'
        load('../../Data/Model Fits/Model Accuracy Fits (Fig 4)/BestFitModel-FourPhases-Older.mat') 
        color = [0,0,128/255];
        titletext = 'Older';
    else 
        error('Improper Input')
    end
    
    %% Plot Spatial Plots of All Four Phases
    figure(400+nn)
    subplot(1,4,1)
    hold on
    for ii = 1:length(yData.LN1.normSubjPx(:,1))
    Cxy = cov(yData.LN1.normSubjPx(ii,:),yData.LN1.normSubjPy(ii,:),1);
    mean_x = mean(yData.LN1.normSubjPx(ii,:));
    mean_y = mean(yData.LN1.normSubjPy(ii,:));
    error_ellipse(Cxy,'mu',[mean_x,mean_y],'conf',0.95,'style','k')
    end
    plot(yData.LN1.PxSubjMean,yData.LN1.PySubjMean,'k','LineWidth',4)
    plot(Young.LN1.Traj.x(:,1),Young.LN1.Traj.x(:,2),'color',color,'LineWidth',4)
    axis equal
    axis([-0.15 0.15 -0.15 0.15])
    title('Late Baseline')

    subplot(1,4,2)
    hold on
    for ii = 1:length(yData.EF1.normSubjPx(:,1))
    Cxy = cov(yData.EF1.normSubjPx(ii,:),yData.EF1.normSubjPy(ii,:),1);
    mean_x = mean(yData.EF1.normSubjPx(ii,:));
    mean_y = mean(yData.EF1.normSubjPy(ii,:));
    error_ellipse(Cxy,'mu',[mean_x,mean_y],'conf',0.95,'style','k')
    end
    plot(yData.EF1.PxSubjMean,yData.EF1.PySubjMean,'k','LineWidth',4)
    plot(Young.EF1.Traj.x(:,1),Young.EF1.Traj.x(:,2),'color',color,'LineWidth',4)
    axis equal
    axis([-0.15 0.15 -0.15 0.15]) 
    title('Early Learning')

    subplot(1,4,3)
    hold on
    for ii = 1:length(yData.LF2.normSubjPx(:,1))
    Cxy = cov(yData.LF2.normSubjPx(ii,:),yData.LF2.normSubjPy(ii,:),1);
    mean_x = mean(yData.LF2.normSubjPx(ii,:));
    mean_y = mean(yData.LF2.normSubjPy(ii,:));
    error_ellipse(Cxy,'mu',[mean_x,mean_y],'conf',0.95,'style','k')
    end
    plot(yData.LF2.PxSubjMean,yData.LF2.PySubjMean,'k','LineWidth',4)
    plot(Young.LF2.Traj.x(:,1),Young.LF2.Traj.x(:,2),'color',color,'LineWidth',4)
    axis equal
    axis([-0.15 0.15 -0.15 0.15])
    title('Late Learning')

    subplot(1,4,4)
    hold on
    for ii = 1:length(yData.EN2.normSubjPx(:,1))
        Cxy = cov(yData.EN2.normSubjPx(ii,:),yData.EN2.normSubjPy(ii,:),1);
        mean_x = mean(yData.EN2.normSubjPx(ii,:));
        mean_y = mean(yData.EN2.normSubjPy(ii,:));
        error_ellipse(Cxy,'mu',[mean_x,mean_y],'conf',0.95,'style','k')
    end
    plot(yData.EN2.PxSubjMean,yData.EN2.PySubjMean,'k','LineWidth',4)
    plot(Young.EN2.Traj.x(:,1),Young.EN2.Traj.x(:,2),'color',color,'LineWidth',4)
    axis equal
    axis([-0.15 0.15 -0.15 0.15])
    title('Early Washout')
    sgtitle(titletext)    
    
    % Toggle ON if wanting to copy / import to illustrator
%     print('-painters','-dsvg',['Fig4-',titletext,'-Spatial'])

    %% Plot Error Metrics Plots of All Four Phases
    phaseNames = {'LN1','EF1','LF2','EN2'};
    phaseNamesPaper = {'LB','EL','LL','EW'};
    legLabel = {'Data','Best'};

    colors = [0.8,0.8,0.8;color];
    
    %% Max Perp Error
    % Of all trials, and standard deviations
    maxPerp = zeros(length(phaseNames),2);
    maxPerpSE = zeros(length(phaseNames),2);
    for jj = 1:length(phaseNames)
        % Of Data
        eval(['ind = find(abs(yData.',phaseNames{jj},...
            '.PxSubjMean) == max(abs(yData.',phaseNames{jj},...
            '.PxSubjMean)),1,''first'');']);
        eval(['maxPerp(jj,1) = yData.',phaseNames{jj},...
            '.PxSubjMean(ind);']);
        eval(['maxPerpSE(jj,1) = 1.96*yData.',phaseNames{jj},...
            '.PxSubjSE(ind);']);

        % Bootstrapped... if using this method
%         eval(['maxPerp(jj,1) = yBootData.',phaseNames{jj},...
%             '.MaxPxMean;']);
%         eval(['maxPerpSE(jj,1) = 1.96*yBootData.',phaseNames{jj},...
%             '.MaxPxSE;']);  

%         % Of Best Fit Model
        eval(['ind = find(abs(Young.',phaseNames{jj},...
            '.Traj.x(:,1)) == max(abs(Young.',...
            phaseNames{jj},'.Traj.x(:,1))),1);']);
        eval(['maxPerp(jj,2) = Young.',phaseNames{jj},'.Traj.x(ind,1);'])
    end

    %Plot Bar
    figure(402+nn)
    sgtitle(titletext)
    subplot(1,3,1)
    hold on
    h = bar(maxPerp,'BarWidth', 1);
    XDATA = bsxfun(@plus, h(1).XData, [h.XOffset]');

    numgroups = size(maxPerp, 1); 
    numbars = size(maxPerp, 2); 
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:numbars
        h(i).FaceColor = colors(i,:);
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        errorbar(x, maxPerp(:,i), maxPerpSE(:,i), 'k', 'linestyle', 'none');
    end
    ylim([-0.1 0.1])
    set(gca,'XTick',1:length(phaseNamesPaper))
    set(gca,'XTickLabel',phaseNamesPaper)
    set(gca,'TickDir','Out')
    ylabel('meters')
    title('Maximum Perpendicular Error')
    legend(h,legLabel)

    %% Estimated Gain 

    % First, Need LN1 Chan and EF1 Chan
    % LN1
    LN1ChanParams = Young.LN1.Params;
    LN1ChanParams.trialLength = length(yDataChan.LN1.PxSubjMean);
    LN1ChanParams.x0 = [yDataChan.LN1.PxSubjMean(1);
                        yDataChan.LN1.PySubjMean(1);
                        yDataChan.LN1.VxSubjMean(1);
                        yDataChan.LN1.VySubjMean(1);
                        yDataChan.LN1.FxSubjMean(1);
                        yDataChan.LN1.FySubjMean(1);
                        yDataChan.LN1.PxSubjMean(end);
                        yDataChan.LN1.PySubjMean(end);0;0];
    SOL_LN1Chan = generateReachTrajectory_Channel(LN1ChanParams);
    ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_LN1Chan.x(:,1)';SOL_LN1Chan.x(:,2)';SOL_LN1Chan.x(:,3)';SOL_LN1Chan.x(:,4)'];
    SOL_LN1Chan.Fx = ForceMat(1,:)';
    SOL_LN1Chan.Fy = ForceMat(2,:)';
    SOL_LN1Chan.AdaptInd= fitgain(SOL_LN1Chan.Fx,SOL_LN1Chan.x(:,4)); 
    Young.LN1Chan.Params = LN1ChanParams;
    Young.LN1Chan.Traj = SOL_LN1Chan;
    %EF1
    EF1ChanParams = Young.EF1.Params;
    EF1ChanParams.trialLength = length(yDataChan.EF1.PxSubjMean);
    EF1ChanParams.x0 = [yDataChan.EF1.PxSubjMean(1);
                        yDataChan.EF1.PySubjMean(1);
                        yDataChan.EF1.VxSubjMean(1);
                        yDataChan.EF1.VySubjMean(1);
                        yDataChan.EF1.FxSubjMean(1);
                        yDataChan.EF1.FySubjMean(1);
                        yDataChan.EF1.PxSubjMean(end);
                        yDataChan.EF1.PySubjMean(end);0;0];
    SOL_EF1Chan = generateReachTrajectory_Channel(EF1ChanParams);
    ForceMat = [-1*stiffConst 0 -1*dampConst 0; 0,0,0,0]*[SOL_EF1Chan.x(:,1)';SOL_EF1Chan.x(:,2)';SOL_EF1Chan.x(:,3)';SOL_EF1Chan.x(:,4)'];
    SOL_EF1Chan.Fx = ForceMat(1,:)';
    SOL_EF1Chan.Fy = ForceMat(2,:)';
    SOL_EF1Chan.AdaptInd= fitgain(SOL_EF1Chan.Fx,SOL_EF1Chan.x(:,4)); 
    Young.EF1Chan.Params = EF1ChanParams;
    Young.EF1Chan.Traj = SOL_EF1Chan;

    fitGain = zeros(length(phaseNames),2);
    fitGainSE = zeros(length(phaseNames),2);
    for jj = 1:length(phaseNames)
        % Of Data
        eval(['fitGain(jj,1) = yDataChan.',phaseNames{jj},'.estGainMean;'])
        eval(['fitGainSE(jj,1) = 1.96*yDataChan.',phaseNames{jj},'.estGainSE;'])
%         % Of Boostrapped Data
%         eval(['fitGain(jj,1) = yBootData.',phaseNames{jj},'.AdaptIndMean;'])
%         eval(['fitGainSE(jj,1) = 1.96*yBootData.',phaseNames{jj},'.AdaptIndSE;'])

        % Of Best Fit Model
        eval(['fitGain(jj,2) = Young.',phaseNames{jj},'Chan.Traj.AdaptInd;'])

    end

    subplot(1,3,2)
    hold on
%     h = bar(-1*fitGain,'BarWidth', 1);
    h = bar(fitGain/-20,'BarWidth', 1);
    XDATA = bsxfun(@plus, h(1).XData, [h.XOffset]');

    numgroups = size(fitGain, 1); 
    numbars = size(fitGain, 2); 
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:numbars
        h(i).FaceColor = colors(i,:);
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
%         errorbar(x, -1*fitGain(:,i), fitGainSE(:,i), 'k', 'linestyle', 'none');
        errorbar(x, fitGain(:,i)./-20, fitGainSE(:,i)./20, 'k', 'linestyle', 'none');
    end
    
    set(gca,'XTick',1:length(phaseNamesPaper))
    set(gca,'XTickLabel',phaseNamesPaper)
    set(gca,'TickDir','Out')
    ylabel('N-s/m')
    ylim([-.25 1])
    title('Curl Gain Estimate')
    legend(h,legLabel)

    %% Max Perp Force
    maxFx = zeros(length(phaseNames),2);
    maxFxSE = zeros(length(phaseNames),2);
    for jj = 1:length(phaseNames)
        % Of Data
        eval(['[~,ind] = max(abs(yDataChan.',phaseNames{jj},'.FxSubjMean));'])
        eval(['maxFx(jj,1) = yDataChan.',phaseNames{jj},'.FxSubjMean(ind);']);
        eval(['maxFxSE(jj,1) = 1.96*yDataChan.',phaseNames{jj},'.FxSubjSE(ind);']);

        % Of Bootstrapped Data
%         eval(['maxFx(jj,1) = yBootData.',phaseNames{jj},'.MaxFxMean;']);
%         eval(['maxFxSE(jj,1) = 1.96*yBootData.',phaseNames{jj},'.MaxFxSE;']);

        % Of Best Fit Model
        eval(['[~,ind] = max(abs(Young.',phaseNames{jj},'Chan.Traj.Fx));']);
        eval(['maxFx(jj,2) = Young.',phaseNames{jj},'Chan.Traj.Fx(ind);'])

    end

    subplot(1,3,3)
    hold on
    h = bar(-1*maxFx,'BarWidth', 1);
    % XDATA=get(get(h,'Children'),'XData');
    XDATA = bsxfun(@plus, h(1).XData, [h.XOffset]');

    numgroups = size(maxFx, 1); 
    numbars = size(maxFx, 2); 
    groupwidth = min(0.8, numbars/(numbars+1.5));
    for i = 1:numbars
        h(i).FaceColor = colors(i,:);
        % Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
        x = (1:numgroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);  % Aligning error bar with individual bar
        errorbar(x, -1*maxFx(:,i), maxFxSE(:,i), 'k', 'linestyle', 'none');
    end

    set(gca,'XTick',1:length(phaseNamesPaper))
    set(gca,'XTickLabel',phaseNamesPaper)
    set(gca,'TickDir','Out')
    ylabel('N')
    ylim([-5 20])
    title('Maximum Horizontal Force')
    legend(h,legLabel)
    
    % Toggle ON if wanting to copy / import to illustrator
%     print('-painters','-dsvg',['Fig4-',titletext,'-ErrorMetrics'])

end
