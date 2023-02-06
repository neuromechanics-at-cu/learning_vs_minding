addpath('../../../Scripts/Supporting Functions')
addpath('../Model-Specific Supporting Functions')

ages = {'Y','O'};
fit_dir = '../../Data/Model Fits/Individual Subject Analysis Fits (Fig 8)';


plot_colors = [0,104,56;46,49,146]./255;
%% Consolidate Data
for ii = 1:length(ages)
    %% Load Data
    if ages{ii} == 'Y'
        load('../../Data/Trajectory Data/Processed Data/Younger-AveragedTrajectoryData-Individual.mat','yData')
    else
        load('../../Data/Trajectory Data/Processed Data/Older-AveragedTrajectoryData-Individual.mat','yData')
    end

    % Process all subjects
    subj_nums = 1:length(yData.LF2.SubjMeanTrLength);
    CR = nan(1,length(yData.LF2.SubjMeanTrLength));
    endcost = nan(1,length(yData.LF2.SubjMeanTrLength));
    learning = nan(1,length(yData.LF2.SubjMeanTrLength));
    objfxn = nan(1,length(yData.LF2.SubjMeanTrLength));

    for nn = 1:length(subj_nums)
        d = dir([fit_dir,'/',num2str(ages{ii}),'-subj',num2str(nn),'*']);
        % Throwout Outliers
%         if ~isempty(d)  &&  ~((nn==10 || nn==13) && ii==1) && ~((nn==2 || nn==13 || nn==14) && ii==2)
        if isempty(d)
            disp([ages{ii},', Subj ',num2str(nn),' not processed'])
        else %if ~isempty(d) 
%             disp(['Age ',ages{ii},' Subj ',num2str(subj_nums(nn)),' good.'])

            load([fit_dir,'/',d(1).name],'LF2Params','AlphaSweep')
%             load([fit_dir,'/',d(1).name],'EN2Params','AlphaSweep')

            if AlphaSweep{1}.ObjFxn > 1e-3 % Original, PxPy Only
%             if AlphaSweep{1}.ObjFxn > 1e-2 % For all Metrics
                disp([ages{ii},', Subj ',num2str(nn),' soln unacceptable'])
            else %if AlphaSweep{1}.ObjFxn > 1e-3
                subj_ind_LF2 = find(yData.LF2.SubjMeanTrLength(nn) == yData.LF2.SubjMeanTrLength(~isnan(yData.LF2.SubjMeanTrLength)),1,'first');
                subj_ind_EN2 = find(yData.EN2.SubjMeanTrLength(nn) == yData.EN2.SubjMeanTrLength(~isnan(yData.EN2.SubjMeanTrLength)),1,'first');
    
                [J,CR(nn)] = calcCostRatios(LF2Params);
%                 [J,CR(nn)] = calcCostRatios(EN2Params);
                endcost(nn) = J.kin.J_LF2_end_POS_norm;
                learning(nn) = AlphaSweep{1}.Params.estGainLF2;

                objfxn(nn) = AlphaSweep{1}.ObjFxn;
    
                % Plot Each Individual
                if 0
                    figure(800+nn+20*ii)
                    subplot(1,2,1)
                    hold on
                    plot(yData.LF2.PxSubjNorm(:,subj_ind_LF2),yData.LF2.PySubjNorm(:,subj_ind_LF2),'k','LineWidth',2)
                    plot(AlphaSweep{1}.LF2.SOL.x(:,1),AlphaSweep{1}.LF2.SOL.x(:,2),'color',plot_colors(ii,:),'LineWidth',2)
                    axis equal
                    title('LF2')
                    xlim([-.10 .10])
                    ylim([-.12 .12])

                    subplot(1,2,2)
                    hold on
                    plot(yData.EN2.PxSubjNorm(:,subj_ind_EN2),yData.EN2.PySubjNorm(:,subj_ind_EN2),'k','LineWidth',2)
                    plot(AlphaSweep{1}.EN2.SOL.x(:,1),AlphaSweep{1}.EN2.SOL.x(:,2),'color',plot_colors(ii,:),'LineWidth',2)
                    axis equal
                    title('EN2')
                    xlim([-.10 .10])
                    ylim([-.12 .12])
                    sgtitle([ages{ii},'-',num2str(nn)])
                end
                % Plot Only Example Subjects for Manuscript
                if 1 && ((nn==4 && ii==1) || (nn==8 && ii==2))
                    figure(800+nn+20*ii)
                    subplot(1,2,1)
                    hold on
                    plot(yData.LF2.PxSubjNorm(:,subj_ind_LF2),yData.LF2.PySubjNorm(:,subj_ind_LF2),'k','LineWidth',2)
                    plot(AlphaSweep{1}.LF2.SOL.x(:,1),AlphaSweep{1}.LF2.SOL.x(:,2),'color',plot_colors(ii,:),'LineWidth',2)
                    axis equal
                    title('LF2')
                    xlim([-.10 .10])
                    ylim([-.14 .14])

                    subplot(1,2,2)
                    hold on
                    plot(yData.EN2.PxSubjNorm(:,subj_ind_EN2),yData.EN2.PySubjNorm(:,subj_ind_EN2),'k','LineWidth',2)
                    plot(AlphaSweep{1}.EN2.SOL.x(:,1),AlphaSweep{1}.EN2.SOL.x(:,2),'color',plot_colors(ii,:),'LineWidth',2)
                    axis equal
                    title('EN2')
                    xlim([-.10 .10])
                    ylim([-.14 .14])
                    sgtitle([ages{ii},'-',num2str(nn)])
                end
            end
        end
    end
    
    %Calculate Outliers
    outlier_id = isoutlier(CR) | isoutlier(learning);
    CR = CR(~outlier_id);
    learning = learning(~outlier_id);
    out_id{ii} = outlier_id;

    % Print Group Summary Data
    if 1
        figure(801)
        hold on
        hl = histogram(learning/-20,0:0.1:1,'Normalization','pdf');
        hl(1).FaceColor = plot_colors(ii,:)/.9;
        x_new = 0:0.01:1.2;
        pd = fitdist(learning'/-20,'normal');
        y_new = normpdf(x_new,pd.mu,pd.sigma);
        plot(x_new,y_new,'color',plot_colors(ii,:),'LineWidth',2)
        xlim([0.2 1.2])
        ylabel('Probability Density')
        xlabel('Proportion Learned')
        title('Learning')
        learn{ii} = learning(~isnan(learning)); 

        figure(802)
        hold on
        hcr = histogram(log(CR),-8:-1,'Normalization','pdf');
        hcr(1).FaceColor = plot_colors(ii,:)/.9;
        x_new = -8:0.01:-2;
        pd = fitdist(log(CR)','normal');
        y_new = normpdf(x_new,pd.mu,pd.sigma);
        plot(x_new,y_new,'color',plot_colors(ii,:),'LineWidth',2)
        xlim([-8 -2])
        ylim([0 0.6])
        title('Cost Ratio')
        ylabel('Probability Density')
        xlabel('log(Kinematic Cost/Effort Cost)')

        %stats
        cr{ii} = log(CR(~isnan(CR))); 
    end
end

%% Stats

[h_learn,p_learn,ci_learn,stats_learn]  = ttest2(learn{1},learn{2})
[h_cr,p_cr,ci_cr,stats_cr]  = ttest2(cr{1},cr{2})
