%% Model Fits Across a Range of Learning


%% Initialize Parameters
phaseNames = {'LF2','EN2'};
ages = ['Y','O'];

for nn = 1:length(ages);
    clearvars -except ages nn phaseNames

    %% Load Data
    if ages(nn) == 'Y'
        load('../../Data/Model Fits/Model Sensitivity Fits (Fig 4)/ModelSensitivity-Younger.mat')
          color = [0,128/255,0];
          light_color = [144/255,238/255,144/255];


    elseif ages(nn) == 'O'
        load('../../Data/Model Fits/Model Sensitivity Fits (Fig 4)/ModelSensitivity-Older.mat')
          color = [0,0,128/255];
          light_color = [190/255,215/255,1];
    else 
        error('Improper Input')
    end

    %% Max Perp Force
    maxFx = zeros(length(phaseNames),length(alpha)+2);
    maxFxSE = zeros(length(phaseNames),length(alpha)+2);

    for jj = 1:length(phaseNames)
        % Of Data
        eval(['[~,ind] = max(abs(yDataChan.',phaseNames{jj},'.FxSubjMean));'])
        eval(['maxFx(jj,1) = yDataChan.',phaseNames{jj},'.FxSubjMean(ind);']);
        eval(['maxFxSE(jj,1) = 1.96*yDataChan.',phaseNames{jj},'.FxSubjSE(ind);']);

        % Of Best Fit Model
        eval(['[~,ind] = max(abs(Young.',phaseNames{jj},'Chan.Traj.Fx));']);
        eval(['maxFx(jj,2) = Young.',phaseNames{jj},'Chan.Traj.Fx(ind);'])

        % Of Alpha Sweep
        for kk = 1:length(alpha)
            eval(['[~,ind] = max(abs(YoungAlphaSweep(kk).',phaseNames{jj},'Chan.Traj.Fx));']);
            eval(['maxFx(jj,kk+2) = YoungAlphaSweep(kk).',phaseNames{jj},'Chan.Traj.Fx(ind);']);
        end
    end

    %% Adaptation Index
    adaptInd = zeros(length(phaseNames),length(alpha)+2);
    adaptIndSE = zeros(length(phaseNames),length(alpha)+2);

    for jj = 1:length(phaseNames)
        % Of Data
        eval(['adaptInd(jj,1) = yDataChan.',phaseNames{jj},'.estGainMean;'])
        eval(['adaptIndSE(jj,1) = 1.96*yDataChan.',phaseNames{jj},'.estGainSE;'])

        % Of Best Fit
        eval(['adaptInd(jj,2) = Young.',phaseNames{jj},'Chan.Traj.AdaptInd;'])

        % Of Alpha Sweep
        for kk = 1:length(alpha)
            eval(['adaptInd(jj,kk+2) = YoungAlphaSweep(kk).',phaseNames{jj},'Chan.Traj.AdaptInd; '])
        end
    end

    %% Max Perp Error
    maxPerp = zeros(length(phaseNames),length(alpha)+2);
    maxPerpSE = zeros(length(phaseNames),length(alpha)+2);
    for jj = 1:length(phaseNames)
        % Of Data
        eval(['[~,ind] = max(abs(yData.',phaseNames{jj},'.PxSubjMean));']);
        eval(['maxPerp(jj,1) = yData.',phaseNames{jj},'.PxSubjMean(ind);']);
        eval(['maxPerpSE(jj,1) = 1.96*yData.',phaseNames{jj},...
            '.PxSubjSE(ind);']);

        % Of Best Fit Model
        eval(['[~,ind] = max(abs(Young.',phaseNames{jj},'.Traj.x(:,1)));'])
        eval(['maxPerp(jj,2) = Young.',phaseNames{jj},'.Traj.x(ind,1);'])

        % Of Alpha Sweep
        for kk = 1:length(alpha)
            eval(['[~,ind] = max(abs(YoungAlphaSweep(kk).',phaseNames{jj},'.Traj.x(:,1)));']);
            eval(['maxPerp(jj,kk+2) = YoungAlphaSweep(kk).',phaseNames{jj},'.Traj.x(ind,1);']);
        end
    end

    %% STATS
    % Make a Binary Matrix of Results
    stats = zeros(length(alpha),3,length(phaseNames));
    for ii = 1:length(alpha)
        for jj = 1:2
            % MaxFx
            if (maxFx(jj,ii+2) < maxFx(jj,1) + maxFxSE(jj,1)) & (maxFx(jj,ii+2) > maxFx(jj,1) - maxFxSE(jj,1))
                stats(ii,1,jj) = 1;
            end
            % AdaptInd
            if (adaptInd(jj,ii+2) < adaptInd(jj,1) + adaptIndSE(jj,1)) & (adaptInd(jj,ii+2) > adaptInd(jj,1) - adaptIndSE(jj,1))
                stats(ii,2,jj) = 1;
            end    
            % MaxPx
            if (maxPerp(jj,ii+2) < maxPerp(jj,1) + maxPerpSE(jj,1)) & (maxPerp(jj,ii+2) > maxPerp(jj,1) - maxPerpSE(jj,1))
                stats(ii,3,jj) = 1;
            end    
        end
    end
    % AND across Phases and Metrics... see what fits all
    stats_byPhase = squeeze(prod(stats,2));
    stats_allPhase_allMetrics = prod(stats_byPhase,2);


    %% MAX PERP FORCE
    figure
    % LF2
    subplot(1,2,1)
    hold on;
    h1t = plot(alpha(logical(stats_allPhase_allMetrics)),-1*maxFx(1,logical([0;0;stats_allPhase_allMetrics])),'o');
    set(h1t,'MarkerEdgeColor','none','MarkerFaceColor',color)
    h1f = plot(alpha(~logical(stats_allPhase_allMetrics)),-1*maxFx(1,logical([0;0;~logical(stats_allPhase_allMetrics)])),'o');
    set(h1f,'MarkerEdgeColor','none','MarkerFaceColor',light_color)    
    bar(0.2,-1*maxFx(1,1),0.1,'FaceColor','k')
    errorbar(0.2,-1*maxFx(1,1),maxFxSE(1,1),'Color',[0.9,0.9,0.9])
    axis([0, 1.2, 0, 20])
    title(phaseNames{1})
    set(gca,'XTick',linspace(0,1.2,13))
    xlab = get(gca,'XTickLabel');
    xlab{1} = '';
    xlab{2} = '';
    xlab{3} = 'Data';
    xlab{4} = '';
    xlab{6} = '';
    xlab{8} = '';
    xlab{10} = '';
    xlab{12} = '';
    set(gca,'XTickLabel',xlab)
    set(gca,'TickDir','out')
    xlabel('LF2 Proportion Learned')
    ylabel('Maximum Horizontal Force (N)')

    % EN2
    subplot(1,2,2)
    hold on;
    h2t = plot(alpha(logical(stats_allPhase_allMetrics)),-1*maxFx(2,logical([0;0;stats_allPhase_allMetrics])),'o');
    set(h2t,'MarkerEdgeColor','none','MarkerFaceColor',color)
    h2f = plot(alpha(~logical(stats_allPhase_allMetrics)),-1*maxFx(2,logical([0;0;~logical(stats_allPhase_allMetrics)])),'o');
    set(h2f,'MarkerEdgeColor','none','MarkerFaceColor',light_color)  
    bar(0.2,-1*maxFx(2,1),0.1,'FaceColor','k')
    errorbar(0.2,-1*maxFx(2,1),maxFxSE(2,1),'Color',[0.9,0.9,0.9])
    axis([0, 1.2, 0, 20])
    title(phaseNames{2})
    set(gca,'XTick',linspace(0,1.2,13))
    xlab = get(gca,'XTickLabel');
    xlab{1} = '';
    xlab{2} = '';
    xlab{3} = 'Data';
    xlab{4} = '';
    xlab{6} = '';
    xlab{8} = '';
    xlab{10} = '';
    xlab{12} = '';
    set(gca,'XTickLabel',xlab)
    set(gca,'TickDir','out')
    xlabel('LF2 Proportion Learned')
    ylabel('Maximum Horizontal Force (N)')

    %% ADAPTATION INDEX
    figure
    % LF2
    subplot(1,2,1)
    hold on
    h1t = plot(alpha(logical(stats_allPhase_allMetrics)),-1*adaptInd(1,logical([0;0;stats_allPhase_allMetrics])),'o');
    set(h1t,'MarkerEdgeColor','none','MarkerFaceColor',color)
    h1f = plot(alpha(~logical(stats_allPhase_allMetrics)),-1*adaptInd(1,logical([0;0;~logical(stats_allPhase_allMetrics)])),'o');
    set(h1f,'MarkerEdgeColor','none','MarkerFaceColor',light_color)  
    bar(0.2,-1*adaptInd(1,1),0.1,'FaceColor','k')
    errorbar(0.2,-1*adaptInd(1,1),adaptIndSE(1,1),'Color',[0.9,0.9,0.9])
    axis([0, 1.2, 0, 24])
    title(phaseNames{1})
    set(gca,'XTick',linspace(0,1.2,13))
    xlab = get(gca,'XTickLabel');
    xlab{1} = '';
    xlab{2} = '';
    xlab{3} = 'Data';
    xlab{4} = '';
    xlab{6} = '';
    xlab{8} = '';
    xlab{10} = '';
    xlab{12} = '';
    set(gca,'XTickLabel',xlab)
    set(gca,'TickDir','out')
    xlabel('LF2 Proportion Learned')
    ylabel('Adaptation Index - Curl Gain Estimate (N-s/m)')

    % EN2
    subplot(1,2,2)
    hold on;
    h2t = plot(alpha(logical(stats_allPhase_allMetrics)),-1*adaptInd(2,logical([0;0;stats_allPhase_allMetrics])),'o');
    set(h2t,'MarkerEdgeColor','none','MarkerFaceColor',color)
    h2f = plot(alpha(~logical(stats_allPhase_allMetrics)),-1*adaptInd(2,logical([0;0;~logical(stats_allPhase_allMetrics)])),'o');
    set(h2f,'MarkerEdgeColor','none','MarkerFaceColor',light_color)  
    bar(0.2,-1*adaptInd(2,1),0.1,'FaceColor','k')
    errorbar(0.2,-1*adaptInd(2,1),adaptIndSE(2,1),'Color',[0.9,0.9,0.9])
    axis([0, 1.2, 0, 24])
    title(phaseNames{2})
    set(gca,'XTick',linspace(0,1.2,13))
    xlab = get(gca,'XTickLabel');
    xlab{1} = '';
    xlab{2} = '';
    xlab{3} = 'Data';
    xlab{4} = '';
    xlab{6} = '';
    xlab{8} = '';
    xlab{10} = '';
    xlab{12} = '';
    set(gca,'XTickLabel',xlab)
    set(gca,'TickDir','out')
    xlabel('LF2 Proportion Learned')
    ylabel('Adaptation Index - Curl Gain Estimate (N-s/m)')

    %% MAX PERP ERROR
    figure
    % LF2
    subplot(1,2,1)
    hold on;
    h1t = plot(alpha(logical(stats_allPhase_allMetrics)),maxPerp(1,logical([0;0;stats_allPhase_allMetrics])),'o');
    set(h1t,'MarkerEdgeColor','none','MarkerFaceColor',color)
    h1f = plot(alpha(~logical(stats_allPhase_allMetrics)),maxPerp(1,logical([0;0;~logical(stats_allPhase_allMetrics)])),'o');
    set(h1f,'MarkerEdgeColor','none','MarkerFaceColor',light_color)  
    bar(0.2,maxPerp(1,1),0.1,'FaceColor','k')
    errorbar(0.2,maxPerp(1,1),maxPerpSE(1,1),'Color',[0.9,0.9,0.9])
    axis([0, 1.2, -0.04, 0.02])
    title(phaseNames{1})
    set(gca,'XTick',linspace(0,1.2,13))
    xlab = get(gca,'XTickLabel');
    xlab{1} = '';
    xlab{2} = '';
    xlab{3} = 'Data';
    xlab{4} = '';
    xlab{6} = '';
    xlab{8} = '';
    xlab{10} = '';
    xlab{12} = '';
    set(gca,'XTickLabel',xlab)
    set(gca,'TickDir','out')
    xlabel('LF2 Proportion Learned')
    ylabel('Maximum Perpendicular Error (m)')

    % EN2
    subplot(1,2,2)
    hold on;
    h2t = plot(alpha(logical(stats_allPhase_allMetrics)),maxPerp(2,logical([0;0;stats_allPhase_allMetrics])),'o');
    set(h2t,'MarkerEdgeColor','none','MarkerFaceColor',color)
    h2f = plot(alpha(~logical(stats_allPhase_allMetrics)),maxPerp(2,logical([0;0;~logical(stats_allPhase_allMetrics)])),'o');
    set(h2f,'MarkerEdgeColor','none','MarkerFaceColor',light_color)  
    bar(0.2,maxPerp(2,1),0.1,'FaceColor','k')
    errorbar(0.2,maxPerp(2,1),maxPerpSE(2,1),'Color',[0.9,0.9,0.9])
    axis([0, 1.2, 0, 0.1])
    title(phaseNames{2})
    set(gca,'XTick',linspace(0,1.2,13))
    xlab = get(gca,'XTickLabel');
    xlab{1} = '';
    xlab{2} = '';
    xlab{3} = 'Data';
    xlab{4} = '';
    xlab{6} = '';
    xlab{8} = '';
    xlab{10} = '';
    xlab{12} = '';
    set(gca,'XTickLabel',xlab)
    set(gca,'TickDir','out')
    xlabel('LF2 Proportion Learned')
    ylabel('Maximum Perpendicular Error (m)')


    %% Spatial Plots
    figure
    subplot(1,2,1)
    hold on
    for ii = 1:length(yData.LF2.normSubjPx(:,1))
        Cxy = cov(yData.LF2.normSubjPx(ii,:),yData.LF2.normSubjPy(ii,:),1);
        mean_x = mean(yData.LF2.normSubjPx(ii,:));
        mean_y = mean(yData.LF2.normSubjPy(ii,:));
        error_ellipse(Cxy,'mu',[mean_x,mean_y],'conf',0.95,'style','k')
    end
    plot(yData.LF2.PxSubjMean,yData.LF2.PySubjMean,'k','LineWidth',4)
    axis equal

    for kk = 1:length(alpha)    
        if ~stats_allPhase_allMetrics(kk)
            plot(YoungAlphaSweep(kk).LF2.Traj.x(:,1),YoungAlphaSweep(kk).LF2.Traj.x(:,2),'Color',light_color,'LineWidth',2)
        end
    end
    for kk = 1:length(alpha)
        if stats_allPhase_allMetrics(kk)
            plot(YoungAlphaSweep(kk).LF2.Traj.x(:,1),YoungAlphaSweep(kk).LF2.Traj.x(:,2),'Color',color,'LineWidth',2)
        end
    end

    axis equal
    axis([-.15 .15 -.15 .15])
    set(gca,'TickDir','out')
    title('LF2')
    xlabel('(m)')
    ylabel('(m)')

    subplot(1,2,2)
    hold on
    for ii = 1:length(yData.EN2.normSubjPx(:,1))
        Cxy = cov(yData.EN2.normSubjPx(ii,:),yData.EN2.normSubjPy(ii,:),1);
        mean_x = mean(yData.EN2.normSubjPx(ii,:));
        mean_y = mean(yData.EN2.normSubjPy(ii,:));
        error_ellipse(Cxy,'mu',[mean_x,mean_y],'conf',0.95,'style','k')
    end
    plot(yData.EN2.PxSubjMean,yData.EN2.PySubjMean,'k','LineWidth',4)
    axis equal


    for kk = 1:length(alpha)    
        if ~stats_allPhase_allMetrics(kk)
            plot(YoungAlphaSweep(kk).EN2.Traj.x(:,1),YoungAlphaSweep(kk).EN2.Traj.x(:,2),'Color',light_color,'LineWidth',2)
        end
    end
    for kk = 1:length(alpha)
        if stats_allPhase_allMetrics(kk)
            plot(YoungAlphaSweep(kk).EN2.Traj.x(:,1),YoungAlphaSweep(kk).EN2.Traj.x(:,2),'Color',color,'LineWidth',2)
        end
    end


    axis equal
    axis([-.15 .15 -.15 .15])
    set(gca,'TickDir','out')
    title('EN2')
    xlabel('(m)')
    ylabel('(m)')
end