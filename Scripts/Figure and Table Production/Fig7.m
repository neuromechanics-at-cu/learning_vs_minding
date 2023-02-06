%% Differences in Subjective Costs


%% Load Data
subjNames = {'Y','O'};
runNames = 1:10;
alphaNames = 0.6:0.05:0.8;
% colors = [1,0,0;1,127/255,0];
colors = [0,128/255,0;0,0,128/255];

load('../../Data/Model Fits/Cost Weight Analysis Fits (Fig 7)/CostRatios_From_Fits.mat')

%% Figure 7: Ratio of Normalized Kinematic Costs to Normalized Effort Costs

figure(701)
% Plot median line
h(1) = semilogy(alphaNames,median(squeeze(err2force(1,:,:))),'color',colors(1,:));
hold on
h(2) = semilogy(alphaNames,median(squeeze(err2force(2,:,:))),'color',colors(2,:));

% Plot Individual points
for ii = 1:length(alphaNames)
    semilogy(alphaNames(ii)*ones(1,length(runNames)),squeeze(err2force(1,:,ii)),'.','MarkerSize',20,'color',colors(1,:))
    semilogy(alphaNames(ii)*ones(1,length(runNames)),squeeze(err2force(2,:,ii)),'.','MarkerSize',20,'color',colors(2,:))
end

% xlim([0.6 0.8])
xlim([min(alphaNames) max(alphaNames)])
xlabel('Proportion Learned')
ylabel('J_{kin,norm}/J_{eff,norm}')
title('Ratio of Normalized Kinematic Costs to Normalized Effort Costs')
legend(h,{'Young','Old'})
set(gca,'XTick',alphaNames);


%% Stats
% Reshape
alphaMat = repmat(alphaNames,10,1);
vecAlpha = reshape(alphaMat,1,50);
vecErr2Force_Y = reshape(squeeze(err2force(1,:,:)),1,50);
vecErr2Force_O = reshape(squeeze(err2force(2,:,:)),1,50);

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


%% Check the Extremes Using Regression
% If cost ratios are equal, is there a solution that describes both younger
% and older adult data? 
CR_Y = log(vecErr2Force_Y)';
alp_Y = vecAlpha';
tbl_Y = table(alp_Y,CR_Y);
mdl_Y = fitlm(tbl_Y,'CR_Y ~ alp_Y');

CR_O = log(vecErr2Force_O)';
alp_O = vecAlpha';
tbl_O = table(alp_O,CR_O);
mdl_O = fitlm(tbl_O,'CR_O ~ alp_O');

int_Y = mdl_Y.Coefficients.Estimate(1);
slope_Y = mdl_Y.Coefficients.Estimate(2);
int_O = mdl_O.Coefficients.Estimate(1);
slope_O = mdl_O.Coefficients.Estimate(2);

% solve for a_Y = 0.9
a_Y_1 = 0.9;
a_O_1 = (slope_Y*a_Y_1 + int_Y - int_O)/slope_O; 

% solve for a_O = 0.6
a_O_2 = 0.6;
a_Y_2 = (slope_O*a_O_2 + int_O - int_Y)/slope_Y;
