%% Process Channel Trial Data from Old or Young Subects
% 
% Throws out bad trials/unnecessary data
% Processes averages for each subject group and phase of the experiment

function D = processRawCurlTrialData(age,trial_bin_type)
% 'subj' is 3-letter acronym, ex: 'FTN'
% 'age' is age group, ex: 'Y' or 'O'

addpath('../../../Scripts/Supporting Functions')

% Check Age Group
if age == 'Y' || age == 'y'
    dirLoc = '../../../Data/Trajectory Data/Raw Data/Younger/RESULTS';
    dirLoc2 = '../../../Data/Trajectory Data/Raw Data/Younger';
    subjs = {'FFO','FFT','FSH','FST','FTN','ITK','OJL','OTG','OUW','SXI',...
            'YEC','ZEH','ZSE'}; 
%     subjs = {'FFO','FFT','FSH','FST','FTN','ITK','OJL','OTG','OUW','SXI',...
%             'YEC','ZEH','ZSE','XZL'}; 
        
elseif age == 'O' || age == 'o'
    dirLoc = '../../../Data/Trajectory Data/Raw Data/Older/RESULTS';
    dirLoc2 = '../../../Data/Trajectory Data/Raw Data/Older';
    subjs = {'BET','EXM','EZC','LGO','MCS','MFQ','MTT','OLT','OSG','OZH',...
            'OZH','QFE','SST','TDS','UTI'};
%     subjs = {'BET','EXM','EZC','LGO','MCS','MFQ','MTT','OLT','OSG','OZH',...
%             'OZT','QFE','SST','TDS','UTI'};
else
    disp('ERROR: Invalid Age Group')
    return
end
    
phaseNames = {'EN1','LN1','EF1','LF1','EF2','LF2','EN2','LN2'};% Defined By Helen

if trial_bin_type == 1
    % First Five
    phaseTrialNums.EN1 = 1:5;
    phaseTrialNums.EF1 = 201:205;
    phaseTrialNums.EF2 = 451:455;
    phaseTrialNums.EN2 = 701:705;
    
    phaseTrialNums.LN1 = 196:200;
    phaseTrialNums.LF1 = 446:450;
    phaseTrialNums.LF2 = 696:700;
    phaseTrialNums.LN2 = 896:900;
elseif trial_bin_type == 2
    % First Ten
    phaseTrialNums.EN1 = 1:10;
    phaseTrialNums.EF1 = 201:210;
    phaseTrialNums.EF2 = 451:460;
    phaseTrialNums.EN2 = 701:710;
    
    phaseTrialNums.LN1 = 191:200;
    phaseTrialNums.LF1 = 441:450;
    phaseTrialNums.LF2 = 691:700;
    phaseTrialNums.LN2 = 891:900;
elseif trial_bin_type == 3
    % Second Five
    phaseTrialNums.EN1 = 6:10;
    phaseTrialNums.EF1 = 206:210;
    phaseTrialNums.EF2 = 456:460;
    phaseTrialNums.EN2 = 706:710;
    
    phaseTrialNums.LN1 = 191:195;
    phaseTrialNums.LF1 = 441:445;
    phaseTrialNums.LF2 = 691:695;
    phaseTrialNums.LN2 = 891:895;
elseif trial_bin_type == 4 
    % 5 for Early, 10 for late
    phaseTrialNums.EN1 = 1:5;
    phaseTrialNums.EF1 = 201:205;
    phaseTrialNums.EF2 = 451:455;
    phaseTrialNums.EN2 = 701:705;
    
    phaseTrialNums.LN1 = 191:200;
    phaseTrialNums.LF1 = 441:450;
    phaseTrialNums.LF2 = 691:700;
    phaseTrialNums.LN2 = 891:900;
elseif trial_bin_type == 5 
    % 5 for Early, 15 for late
    phaseTrialNums.EN1 = 1:5;
    phaseTrialNums.EF1 = 201:205;
    phaseTrialNums.EF2 = 451:455;
    phaseTrialNums.EN2 = 701:705;
    
    phaseTrialNums.LN1 = 186:200;
    phaseTrialNums.LF1 = 436:450;
    phaseTrialNums.LF2 = 686:700;
    phaseTrialNums.LN2 = 886:900;
else
    error('incorrect trial bin type')
end

for k = 1:length(subjs)
    clear TR
    % Check/Load Subject Data
    if age == 'Y' || age == 'y'
        if exist([dirLoc2,'/',subjs{k},'.mat'])
            load([dirLoc2,'/',subjs{k},'.mat'],'T')
        else
            disp('ERROR: Invalid Subject')
            return
        end
    elseif age == 'O' || age == 'o'
        if exist([dirLoc2,'/',subjs{k},'.mat'])
            load([dirLoc2,'/',subjs{k},'.mat'],'T')
        else
            disp('ERROR: Invalid Subject')
            return
        end
    end

    % Trial Numbers/Types
    allTrialNums = find(strcmp(T.trialtypename,'null') | strcmp(T.trialtypename,'curl'));
    allCTrialNums = find(strcmp(T.trialtypename,'clamp'));
    % Eliminate "Pull" Trajectories
    allTrialNums = allTrialNums(rem(allTrialNums,2) == 1);
    allCTrialNums = allCTrialNums(rem(allCTrialNums,2) == 1);
    % Eliminate Bad Trials
    badTrialNums = [];
    goodTrialNums = [];
    for j = 1:900
%         if T.framedata(j).statenumber(end) < 7 | T.framedata(j).statenumber(end) > 8
        if T.framedata(j).statenumber(end) < 7
            badTrialNums = [badTrialNums, j];
        else
            goodTrialNums = [goodTrialNums, j];
        end
    end
    
    % Get Trial Numbers for Specified Phase
    for j = 1:length(phaseNames)
%         eval(['trialNums.',phaseNames{j},' = intersect(allTrialNums,phaseTrialNums.',phaseNames{j},');'])
%         eval(['ctrialNums.',phaseNames{j},' = intersect(allCTrialNums,phaseTrialNums.',phaseNames{j},');'])
        eval(['trialNums.',phaseNames{j},' = intersect(intersect(allTrialNums,phaseTrialNums.',phaseNames{j},'),goodTrialNums);'])
        eval(['ctrialNums.',phaseNames{j},' = intersect(intersect(allCTrialNums,phaseTrialNums.',phaseNames{j},'),goodTrialNums);'])

        % Initialize Velocity Filter Trial Structure
        eval([phaseNames{j},'.trialNums{k}.tooLong = [];'])
        eval([phaseNames{j},'.trialNums{k}.tooShort = [];'])
        eval([phaseNames{j},'.trialNums{k}.justRight = [];'])

        % Find Number of Trials
        eval(['nTrials = length(trialNums.',phaseNames{j},');'])
        
        % Find Max Trial Length
        maxLength = 0;
        for m = 1:nTrials
            eval(['lengths(m) = length(T.framedata(1,trialNums.',phaseNames{j},'(m)).x);'])
            if lengths(m) > maxLength
               maxLength = lengths(m); 
            end
        end
        
        % Preallocate Matrix
        eval(['TR.trials.Px.',phaseNames{j},'.out = nan(maxLength,nTrials);'])
        eval(['TR.trials.Py.',phaseNames{j},'.out = nan(maxLength,nTrials);'])
        eval(['TR.trials.Vx.',phaseNames{j},'.out = nan(maxLength,nTrials);'])
        eval(['TR.trials.Vy.',phaseNames{j},'.out = nan(maxLength,nTrials);'])
                
        % Fill Matrix with Trial Data
        for m = 1:nTrials
            eval(['TR.trials.Px.',phaseNames{j},'.out(1:lengths(m),m) = T.framedata(1,trialNums.',phaseNames{j},'(m)).x;'])
            eval(['TR.trials.Py.',phaseNames{j},'.out(1:lengths(m),m) = T.framedata(1,trialNums.',phaseNames{j},'(m)).y;'])
            eval(['TR.trials.Vx.',phaseNames{j},'.out(1:lengths(m),m) = T.framedata(1,trialNums.',phaseNames{j},'(m)).vx;'])
            eval(['TR.trials.Vy.',phaseNames{j},'.out(1:lengths(m),m) = T.framedata(1,trialNums.',phaseNames{j},'(m)).vy;'])
        end
            
        % Add Ax, Ay
        eval(['TR.trials.Ax.',phaseNames{j},'.out = nan(maxLength,nTrials);'])
        eval(['TR.trials.Ay.',phaseNames{j},'.out = nan(maxLength,nTrials);'])
        for m = 1:nTrials
            eval(['TR.trials.Ax.',phaseNames{j},'.out(1:lengths(m),m) = NumDiff_3pt(T.framedata(1,trialNums.',phaseNames{j},'(m)).vx,[1/200:1/200:lengths(m)/200]'');'])
            eval(['TR.trials.Ay.',phaseNames{j},'.out(1:lengths(m),m) = NumDiff_3pt(T.framedata(1,trialNums.',phaseNames{j},'(m)).vy,[1/200:1/200:lengths(m)/200]'');'])
            % Add butterworth, zero shift filter, to smooth acceleration
            [bb,aa] = butter(9,.1,'low'); % 0.2 ~ 20 Hz
            eval(['TR.trials.Ax.',phaseNames{j},'.out(1:lengths(m),m) = filtfilt(bb,aa,TR.trials.Ax.',phaseNames{j},'.out(1:lengths(m),m));'])
            eval(['TR.trials.Ay.',phaseNames{j},'.out(1:lengths(m),m) = filtfilt(bb,aa,TR.trials.Ay.',phaseNames{j},'.out(1:lengths(m),m));'])
        end
            
        
        if nTrials >= 1
        maxLength = 0;
        % Find Movement Start (Vy > 0.03 m/s) and Movement End (Vy < 0.03 m/s)   
        for m = 1:nTrials 
            eval(['start = find(TR.trials.Vy.',phaseNames{j},'.out(:,m) > 0.03,1);'])
            % CMH Changed 5/18/22
%             eval(['start = find(TR.trials.Vy.',phaseNames{j},...
%                 '.out(:,m) > 0.03 & (abs(TR.trials.Px.',phaseNames{j},...
%                 '.out(:,m)) > 0.008 | abs(TR.trials.Py.',phaseNames{j},...
%                 '.out(:,m) + 0.1) > 0.008),1);'])
%             if start>15
%                 start = start - 15;
%             else 
%                 start = 1;
%             end
            eval(['finish = start - 2 + ',...
                'find(TR.trials.Vy.',phaseNames{j},'.out(start:end,m) < 0.03',...
                '& abs(TR.trials.Px.',phaseNames{j},'.out(start:end,m)) < 0.008',...
                '& abs(TR.trials.Py.',phaseNames{j},'.out(start:end,m) - 0.1) < 0.008,1);'])
            % CMH Added 5/22
            eval(['vel_finish = start - 2 + find((TR.trials.Vy.',phaseNames{j},'.out(start:end,m) < 0.03)',...
                                                '& (abs(TR.trials.Vx.',phaseNames{j},'.out(start:end,m)) < 0.03)',...
                                                '& (TR.trials.Py.',phaseNames{j},'.out(start:end,m) > 0.07),1);'])
            if isempty(finish) | ((finish - vel_finish) >= 50)
                finish = vel_finish;
            end

            if isempty(finish)
                eval(['finish = length(TR.trials.Vy.',phaseNames{j},'.out(:,m));'])
            end
            % Resave Trajectories to New Length
            eval(['TR.trials.Px.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.trials.Px.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.trials.Py.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.trials.Py.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.trials.Vx.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.trials.Vx.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.trials.Vy.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.trials.Vy.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.trials.Ax.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.trials.Ax.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.trials.Ay.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.trials.Ay.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.trials.Px.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
            eval(['TR.trials.Py.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
            eval(['TR.trials.Vx.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
            eval(['TR.trials.Vy.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
            eval(['TR.trials.Ax.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
            eval(['TR.trials.Ay.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
        end
        % Remove NaNs
        eval(['[inan,jnan] = find(~isnan(TR.trials.Px.',phaseNames{j},'.out));'])
        maxLength = max(inan);
        eval(['TR.trials.Px.',phaseNames{j},'.out = TR.trials.Px.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.trials.Py.',phaseNames{j},'.out = TR.trials.Py.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.trials.Vx.',phaseNames{j},'.out = TR.trials.Vx.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.trials.Vy.',phaseNames{j},'.out = TR.trials.Vy.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.trials.Ax.',phaseNames{j},'.out = TR.trials.Ax.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.trials.Ay.',phaseNames{j},'.out = TR.trials.Ay.',phaseNames{j},'.out(1:maxLength,:);'])

        end
        
        % Eliminate Botched Trials 
        
        % One from 'BET', EN2, Trial 1, has a "blip" from the robot error
%         if subjs{k} == 'BET' & phaseNames{j} == 'EN2' 
%              eval(['trialNums.',phaseNames{j},'(1) = [];'])
%              eval(['TR.trials.Px.',phaseNames{j},'.out(:,1) = [];'])
%              eval(['TR.trials.Py.',phaseNames{j},'.out(:,1) = [];'])
%              eval(['TR.trials.Vx.',phaseNames{j},'.out(:,1) = [];'])
%              eval(['TR.trials.Vy.',phaseNames{j},'.out(:,1) = [];'])
%              eval(['TR.trials.Ax.',phaseNames{j},'.out(:,1) = [];'])
%              eval(['TR.trials.Ay.',phaseNames{j},'.out(:,1) = [];'])
%         end
        
        % (Ex: Subj-'FTN', Phase-EF2)
         minLength = 50;%50     60  (300 ms)
         maxLength = 300;%300   120 (600 ms)
         if eval(['size(TR.trials.Px.',phaseNames{j},'.out,2)==1'])
             if eval(['length(TR.trials.Px.',phaseNames{j},'.out) < minLength'])
                 % CMHNEW
                 eval([phaseNames{j},'.trialNums{k}.tooShort = [',phaseNames{j},'.trialNums{k}.tooShort,trialNums.',phaseNames{j},'];'])
                 %
                 eval(['trialNums.',phaseNames{j},' = [];']) 
                 eval(['TR.trials.Px.',phaseNames{j},'.out = [];'])
                 eval(['TR.trials.Py.',phaseNames{j},'.out = [];'])
                 eval(['TR.trials.Vx.',phaseNames{j},'.out = [];'])
                 eval(['TR.trials.Vy.',phaseNames{j},'.out = [];'])
                 eval(['TR.trials.Ax.',phaseNames{j},'.out = [];'])
                 eval(['TR.trials.Ay.',phaseNames{j},'.out = [];'])
             elseif eval(['length(TR.trials.Px.',phaseNames{j},'.out) > maxLength'])
                 % CMHNEW
                 eval([phaseNames{j},'.trialNums{k}.tooLong = [',phaseNames{j},'.trialNums{k}.tooLong,trialNums.',phaseNames{j},'];'])
                 %
                 eval(['trialNums.',phaseNames{j},' = [];']) 
                 eval(['TR.trials.Px.',phaseNames{j},'.out = [];'])
                 eval(['TR.trials.Py.',phaseNames{j},'.out = [];'])
                 eval(['TR.trials.Vx.',phaseNames{j},'.out = [];'])
                 eval(['TR.trials.Vy.',phaseNames{j},'.out = [];'])
                 eval(['TR.trials.Ax.',phaseNames{j},'.out = [];'])
                 eval(['TR.trials.Ay.',phaseNames{j},'.out = [];'])
             end
             % CMHNEW
             eval([phaseNames{j},'.trialNums{k}.justRight = [',phaseNames{j},'.trialNums{k}.justRight,trialNums.',phaseNames{j},'];'])
             %
         else
             % Remove if Too Short
             eval(['[II,JJ] = find(isnan(TR.trials.Px.',phaseNames{j},'.out));'])
             III = find(II < minLength);
             % CMHNEW
             eval([phaseNames{j},'.trialNums{k}.tooShort = [',phaseNames{j},'.trialNums{k}.tooShort,unique(JJ(III))];'])
             %
             eval(['trialNums.',phaseNames{j},'(JJ(III)) = [];']) 
             eval(['TR.trials.Px.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             eval(['TR.trials.Py.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             eval(['TR.trials.Vx.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             eval(['TR.trials.Vy.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             eval(['TR.trials.Ax.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             eval(['TR.trials.Ay.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             % Remove if Too Long
             eval(['[nII,nJJ] = find(~isnan(TR.trials.Px.',phaseNames{j},'.out));'])
             nIII = find(nII > maxLength);
             % CMHNEW
             eval([phaseNames{j},'.trialNums{k}.tooLong = [',phaseNames{j},'.trialNums{k}.tooLong,unique(nJJ(nIII))];'])
             %
             eval(['trialNums.',phaseNames{j},'(nJJ(nIII)) = [];'])
             eval(['TR.trials.Px.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             eval(['TR.trials.Py.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             eval(['TR.trials.Vx.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             eval(['TR.trials.Vy.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             eval(['TR.trials.Ax.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             eval(['TR.trials.Ay.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             % CMHNEW
             eval([phaseNames{j},'.trialNums{k}.justRight = [',phaseNames{j},'.trialNums{k}.justRight,trialNums.',phaseNames{j},'];'])
             %
         end
         
        % Re-remove NaNs
        eval(['[inan,jnan] = find(~isnan(TR.trials.Px.',phaseNames{j},'.out));'])
        maxLength = max(inan);
        eval(['TR.trials.Px.',phaseNames{j},'.out = TR.trials.Px.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.trials.Py.',phaseNames{j},'.out = TR.trials.Py.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.trials.Vx.',phaseNames{j},'.out = TR.trials.Vx.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.trials.Vy.',phaseNames{j},'.out = TR.trials.Vy.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.trials.Ax.',phaseNames{j},'.out = TR.trials.Ax.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.trials.Ay.',phaseNames{j},'.out = TR.trials.Ay.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.trials.Px.',phaseNames{j},'.mean = [];'])
        eval(['TR.trials.Py.',phaseNames{j},'.mean = [];'])
        eval(['TR.trials.Vx.',phaseNames{j},'.mean = [];'])
        eval(['TR.trials.Vy.',phaseNames{j},'.mean = [];'])
        eval(['TR.trials.Ax.',phaseNames{j},'.mean = [];'])
        eval(['TR.trials.Ay.',phaseNames{j},'.mean = [];'])
        eval(['TR.trials.Px.',phaseNames{j},'.std = [];'])
        eval(['TR.trials.Py.',phaseNames{j},'.std = [];'])
        eval(['TR.trials.Vx.',phaseNames{j},'.std = [];'])
        eval(['TR.trials.Vy.',phaseNames{j},'.std = [];'])
        eval(['TR.trials.Ax.',phaseNames{j},'.std = [];'])
        eval(['TR.trials.Ay.',phaseNames{j},'.std = [];']) 
        % Find Within Subject Mean, STD
        for p = 1:maxLength 
            eval(['TR.trials.Px.',phaseNames{j},'.mean(p,1) = nanmean(TR.trials.Px.',phaseNames{j},'.out(p,:));'])
            eval(['TR.trials.Py.',phaseNames{j},'.mean(p,1) = nanmean(TR.trials.Py.',phaseNames{j},'.out(p,:));'])
            eval(['TR.trials.Vx.',phaseNames{j},'.mean(p,1) = nanmean(TR.trials.Vx.',phaseNames{j},'.out(p,:));'])
            eval(['TR.trials.Vy.',phaseNames{j},'.mean(p,1) = nanmean(TR.trials.Vy.',phaseNames{j},'.out(p,:));'])
            eval(['TR.trials.Ax.',phaseNames{j},'.mean(p,1) = nanmean(TR.trials.Ax.',phaseNames{j},'.out(p,:));'])
            eval(['TR.trials.Ay.',phaseNames{j},'.mean(p,1) = nanmean(TR.trials.Ay.',phaseNames{j},'.out(p,:));'])
            eval(['TR.trials.Px.',phaseNames{j},'.std(p,1) = nanstd(TR.trials.Px.',phaseNames{j},'.out(p,:));'])
            eval(['TR.trials.Py.',phaseNames{j},'.std(p,1) = nanstd(TR.trials.Py.',phaseNames{j},'.out(p,:));'])
            eval(['TR.trials.Vx.',phaseNames{j},'.std(p,1) = nanstd(TR.trials.Vx.',phaseNames{j},'.out(p,:));'])
            eval(['TR.trials.Vy.',phaseNames{j},'.std(p,1) = nanstd(TR.trials.Vy.',phaseNames{j},'.out(p,:));'])
            eval(['TR.trials.Ax.',phaseNames{j},'.std(p,1) = nanstd(TR.trials.Ax.',phaseNames{j},'.out(p,:));'])
            eval(['TR.trials.Ay.',phaseNames{j},'.std(p,1) = nanstd(TR.trials.Ay.',phaseNames{j},'.out(p,:));'])
        end 
        % Re-Find Number of Trials
        eval(['nTrials = length(trialNums.',phaseNames{j},');'])
        
        % Track Max Subject Trial Length
        eval([phaseNames{j},'.SubjTrLength(k) = length(TR.trials.Px.',phaseNames{j},'.mean);'])
        % Track Mean Subject Trial Length
        eval(['trLength = ones(1,nTrials)*size(TR.trials.Px.',phaseNames{j},'.out,1);'])
        for jj = 1:nTrials
            eval(['notLongest = ~isempty(find(isnan(TR.trials.Px.',phaseNames{j},'.out(:,jj)),1,''first''));'])
            if notLongest    
                eval(['trLength(jj) = find(isnan(TR.trials.Px.',phaseNames{j},'.out(:,jj)),1,''first'')-1;'])
            end
        end
        % Not going to use a floor call yet.
        eval([phaseNames{j},'.SubjMeanTrLength(k) = mean(trLength);'])
        
        % Re-sample individual trials to subject mean length, then average
        % for subject mean-normed trajectory
        % Resample Subject Trajectories so the same time length
        avgTrLengthSubj = floor(mean(trLength));
        if isnan(avgTrLengthSubj)
            avgTrLengthSubj = 0;
        end
        PxSubjNorm = zeros(avgTrLengthSubj,nTrials);
        PySubjNorm = zeros(avgTrLengthSubj,nTrials);
        VxSubjNorm = zeros(avgTrLengthSubj,nTrials);
        VySubjNorm = zeros(avgTrLengthSubj,nTrials);
        AxSubjNorm = zeros(avgTrLengthSubj,nTrials);
        AySubjNorm = zeros(avgTrLengthSubj,nTrials);
        for jj = 1:nTrials;
            eval(['PxSubjNorm(:,jj) = resample(TR.trials.Px.',phaseNames{j},...
                '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])
            eval(['PySubjNorm(:,jj) = resample(TR.trials.Py.',phaseNames{j},...
                '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])            
            eval(['VxSubjNorm(:,jj) = resample(TR.trials.Vx.',phaseNames{j},...
                '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])            
            eval(['VySubjNorm(:,jj) = resample(TR.trials.Vy.',phaseNames{j},...
                '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])
            eval(['AxSubjNorm(:,jj) = resample(TR.trials.Ax.',phaseNames{j},...
                '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])            
            eval(['AySubjNorm(:,jj) = resample(TR.trials.Ay.',phaseNames{j},...
                '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])
        end
        eval(['TR.trials.Px.',phaseNames{j},'.normed = [];'])
        eval(['TR.trials.Py.',phaseNames{j},'.normed = [];'])
        eval(['TR.trials.Vx.',phaseNames{j},'.normed = [];'])
        eval(['TR.trials.Vy.',phaseNames{j},'.normed = [];'])
        eval(['TR.trials.Ax.',phaseNames{j},'.normed = [];'])
        eval(['TR.trials.Ay.',phaseNames{j},'.normed = [];'])
        for ii = 1:avgTrLengthSubj
            eval(['TR.trials.Px.',phaseNames{j},'.normed(ii,:) = mean(PxSubjNorm(ii,:));'])
            eval(['TR.trials.Py.',phaseNames{j},'.normed(ii,:) = mean(PySubjNorm(ii,:));'])
            eval(['TR.trials.Vx.',phaseNames{j},'.normed(ii,:) = mean(VxSubjNorm(ii,:));'])
            eval(['TR.trials.Vy.',phaseNames{j},'.normed(ii,:) = mean(VySubjNorm(ii,:));'])
            eval(['TR.trials.Ax.',phaseNames{j},'.normed(ii,:) = mean(AxSubjNorm(ii,:));'])
            eval(['TR.trials.Ay.',phaseNames{j},'.normed(ii,:) = mean(AySubjNorm(ii,:));'])
        end

        % Assemble Subject Mean/STD into Structure
        if k == 1
            eval([phaseNames{j},'.PxSubj = TR.trials.Px.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.PySubj = TR.trials.Py.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.VxSubj = TR.trials.Vx.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.VySubj = TR.trials.Vy.',phaseNames{j},'.mean;']) 
            eval([phaseNames{j},'.AxSubj = TR.trials.Ax.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.AySubj = TR.trials.Ay.',phaseNames{j},'.mean;']) 
        elseif nTrials > 0
            eval(['trialLength = size(TR.trials.Px.',phaseNames{j},'.mean,1);'])
            eval(['matrixLength = size(',phaseNames{j},'.PxSubj,1);'])
            if trialLength > matrixLength
                eval([phaseNames{j},'.PxSubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.PxSubj,2));'])
                eval([phaseNames{j},'.PySubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.PySubj,2));'])
                eval([phaseNames{j},'.VxSubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.VxSubj,2));'])
                eval([phaseNames{j},'.VySubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.VySubj,2));'])  
                eval([phaseNames{j},'.AxSubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.AxSubj,2));'])
                eval([phaseNames{j},'.AySubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.AySubj,2));'])                  
            elseif trialLength < matrixLength
                eval(['TR.trials.Px.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.trials.Py.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.trials.Vx.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.trials.Vy.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.trials.Ax.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.trials.Ay.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
            end
            eval([phaseNames{j},'.PxSubj(:,end+1) = TR.trials.Px.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.PySubj(:,end+1) = TR.trials.Py.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.VxSubj(:,end+1) = TR.trials.Vx.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.VySubj(:,end+1) = TR.trials.Vy.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.AxSubj(:,end+1) = TR.trials.Ax.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.AySubj(:,end+1) = TR.trials.Ay.',phaseNames{j},'.mean;'])
        end
        
        % Assemble Subject NORMED Mean/STD into Structure
        if k == 1
            eval([phaseNames{j},'.PxSubjNorm = TR.trials.Px.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.PySubjNorm = TR.trials.Py.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.VxSubjNorm = TR.trials.Vx.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.VySubjNorm = TR.trials.Vy.',phaseNames{j},'.normed;']) 
            eval([phaseNames{j},'.AxSubjNorm = TR.trials.Ax.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.AySubjNorm = TR.trials.Ay.',phaseNames{j},'.normed;']) 
        elseif nTrials > 0
            eval(['trialLength = size(TR.trials.Px.',phaseNames{j},'.normed,1);'])
            eval(['matrixLength = size(',phaseNames{j},'.PxSubjNorm,1);'])
            if trialLength > matrixLength
                eval([phaseNames{j},'.PxSubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.PxSubjNorm,2));'])
                eval([phaseNames{j},'.PySubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.PySubjNorm,2));'])
                eval([phaseNames{j},'.VxSubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.VxSubjNorm,2));'])
                eval([phaseNames{j},'.VySubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.VySubjNorm,2));'])  
                eval([phaseNames{j},'.AxSubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.AxSubjNorm,2));'])
                eval([phaseNames{j},'.AySubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.AySubjNorm,2));'])                  
            elseif trialLength < matrixLength
                eval(['TR.trials.Px.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.trials.Py.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.trials.Vx.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.trials.Vy.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.trials.Ax.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.trials.Ay.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
            end
            eval([phaseNames{j},'.PxSubjNorm(:,end+1) = TR.trials.Px.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.PySubjNorm(:,end+1) = TR.trials.Py.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.VxSubjNorm(:,end+1) = TR.trials.Vx.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.VySubjNorm(:,end+1) = TR.trials.Vy.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.AxSubjNorm(:,end+1) = TR.trials.Ax.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.AySubjNorm(:,end+1) = TR.trials.Ay.',phaseNames{j},'.normed;'])
        end
        
        if k == 1
            eval([phaseNames{j},'.Px = TR.trials.Px.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Py = TR.trials.Py.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Vx = TR.trials.Vx.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Vy = TR.trials.Vy.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Ax = TR.trials.Ax.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Ay = TR.trials.Ay.',phaseNames{j},'.out;'])
        elseif nTrials > 0
            eval(['trialLength = size(TR.trials.Px.',phaseNames{j},'.out,1);'])
            eval(['matrixLength = size(',phaseNames{j},'.Px,1);'])
            if trialLength > matrixLength
                eval([phaseNames{j},'.Px(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Px,2));'])
                eval([phaseNames{j},'.Py(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Py,2));'])
                eval([phaseNames{j},'.Vx(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Vx,2));'])
                eval([phaseNames{j},'.Vy(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Vy,2));'])            
                eval([phaseNames{j},'.Ax(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Ax,2));'])
                eval([phaseNames{j},'.Ay(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Ay,2));'])            
            elseif trialLength < matrixLength
                eval(['TR.trials.Px.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
                eval(['TR.trials.Py.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
                eval(['TR.trials.Vx.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
                eval(['TR.trials.Vy.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
                eval(['TR.trials.Ax.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
                eval(['TR.trials.Ay.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
            end
            eval([phaseNames{j},'.Px(:,end+1:end+nTrials) = TR.trials.Px.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Py(:,end+1:end+nTrials) = TR.trials.Py.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Vx(:,end+1:end+nTrials) = TR.trials.Vx.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Vy(:,end+1:end+nTrials) = TR.trials.Vy.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Ax(:,end+1:end+nTrials) = TR.trials.Ax.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Ay(:,end+1:end+nTrials) = TR.trials.Ay.',phaseNames{j},'.out;'])            
        end
        
        % Get Curl Field Estimation
        %      Calculates average if phase has multiple catch trials
        estGain = 0;
        eval(['nTemp = length(ctrialNums.',phaseNames{j},');']);
        for i = 1:nTemp
            eval(['estGain = estGain + fitgain(T.framedata(1,ctrialNums.',phaseNames{j},'(i)).fx,T.framedata(1,ctrialNums.',phaseNames{j},'(i)).vy);'])
        end
%         eval([phaseNames{j},'.estGain(k) = -1*estGain/length(ctrialNums.',phaseNames{j},');'])%Takes averages, converts to negative (from abs. value)
        eval([phaseNames{j},'.estGain(k) = estGain/length(ctrialNums.',phaseNames{j},');'])%Takes averages,
                    
        %Max Perp and Angular Error
        if k == 1
           eval([phaseNames{j},'.maxPerp = [];'])
            if nTrials > 0
                for i = 1:nTrials
                    eval(['[~,maxPerpInd] = max(abs(TR.trials.Px.',phaseNames{j},'.out(:,i)));'])
                    eval([phaseNames{j},'.maxPerp = TR.trials.Px.',phaseNames{j},'.out(maxPerpInd,i);'])
                    
                    eval(['[~,maxVelInd] = max(TR.trials.Vx.',phaseNames{j},...
                        '.out(:,i).^2 + TR.trials.Vy.',phaseNames{j},'.out(:,i).^2);'])
                    eval([phaseNames{j},'.angErr = acosd(dot([TR.trials.Px.',phaseNames{j},...
                        '.out(maxVelInd,i) - 0; TR.trials.Py.',phaseNames{j},'.out(maxVelInd,i) + 0.1],[0;1])/(norm([TR.trials.Px.',phaseNames{j},...
                        '.out(maxVelInd,i) - 0; TR.trials.Py.',phaseNames{j},'.out(maxVelInd,i) + 0.1])*norm([0;1])));'])
                end
            end
        else
            eval(['totTrials = length(',phaseNames{j},'.maxPerp);'])
            for i = 1:nTrials
                eval(['[~,maxPerpInd] = max(abs(TR.trials.Px.',phaseNames{j},'.out(:,i)));'])
                eval([phaseNames{j},'.maxPerp(totTrials+i) = TR.trials.Px.',phaseNames{j},'.out(maxPerpInd,i);'])
                
                eval(['[~,maxVelInd] = max(TR.trials.Vx.',phaseNames{j},...
                        '.out(:,i).^2 + TR.trials.Vy.',phaseNames{j},'.out(:,i).^2);'])
                eval([phaseNames{j},'.angErr(totTrials+i) = acosd(dot([TR.trials.Px.',phaseNames{j},...
                        '.out(maxVelInd,i) - 0; TR.trials.Py.',phaseNames{j},'.out(maxVelInd,i) + 0.1],[0;1])/(norm([TR.trials.Px.',phaseNames{j},...
                        '.out(maxVelInd,i) - 0; TR.trials.Py.',phaseNames{j},'.out(maxVelInd,i) + 0.1])*norm([0;1])));'])

            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%
% Average Subject Data

% if length(subjs) > 1
    for j = 1:length(phaseNames)
        % Find Valid Trials for Given Phase
%         eval(['validSubjs = find(~isnan(',phaseNames{j},'.Px(:,1)));'])
%        eval(['validTrials = find(~isnan(',phaseNames{j},'.Px(:,1)));'])
        eval(['numTrials = size(',phaseNames{j},'.Px,2);'])
        
        if numTrials > 0 
        % Find Average Trial Length
        eval(['trLength = ones(1,numTrials)*size(',phaseNames{j},'.Px,1);'])
        for jj = 1:numTrials
            eval(['notLongest = ~isempty(find(isnan(',phaseNames{j},'.Px(:,jj)),1,''first''));'])
            if notLongest    
                eval(['trLength(jj) = find(isnan(',phaseNames{j},'.Px(:,jj)),1,''first'')-1;'])
            end
        end
        % Is Floor a bad assumption? Would Rounding up be better or worse?
        avgTrLength = floor(mean(trLength));

        % Resample Trajectories so the same time length
        normPx = zeros(avgTrLength,numTrials);
        normPy = zeros(avgTrLength,numTrials);
        normVx = zeros(avgTrLength,numTrials);
        normVy = zeros(avgTrLength,numTrials);
        normAx = zeros(avgTrLength,numTrials);
        normAy = zeros(avgTrLength,numTrials);
        for jj = 1:numTrials;
            eval(['normPx(:,jj) = resample(',phaseNames{j},...
                '.Px(1:trLength(jj),jj),avgTrLength,trLength(jj),0);'])
            eval(['normPy(:,jj) = resample(',phaseNames{j},...
                '.Py(1:trLength(jj),jj),avgTrLength,trLength(jj),0);'])            
            eval(['normVx(:,jj) = resample(',phaseNames{j},...
                '.Vx(1:trLength(jj),jj),avgTrLength,trLength(jj),0);'])            
            eval(['normVy(:,jj) = resample(',phaseNames{j},...
                '.Vy(1:trLength(jj),jj),avgTrLength,trLength(jj),0);'])
            eval(['normAx(:,jj) = resample(',phaseNames{j},...
                '.Ax(1:trLength(jj),jj),avgTrLength,trLength(jj),0);'])            
            eval(['normAy(:,jj) = resample(',phaseNames{j},...
                '.Ay(1:trLength(jj),jj),avgTrLength,trLength(jj),0);'])
        end
    
                
        %Save Normalized Trajectories
        eval([phaseNames{j},'.normPx = normPx;'])
        eval([phaseNames{j},'.normPy = normPy;'])
        eval([phaseNames{j},'.normVx = normVx;'])
        eval([phaseNames{j},'.normVy = normVy;'])
        eval([phaseNames{j},'.normAx = normAx;'])
        eval([phaseNames{j},'.normAy = normAy;'])
        
        % Find the average trajectory
        for ii = 1:avgTrLength
            eval([phaseNames{j},'.PxMean(ii) = mean(normPx(ii,:));'])
            eval([phaseNames{j},'.PyMean(ii) = mean(normPy(ii,:));'])
            eval([phaseNames{j},'.VxMean(ii) = mean(normVx(ii,:));'])
            eval([phaseNames{j},'.VyMean(ii) = mean(normVy(ii,:));'])
            eval([phaseNames{j},'.AxMean(ii) = mean(normAx(ii,:));'])
            eval([phaseNames{j},'.AyMean(ii) = mean(normAy(ii,:));'])            
            eval([phaseNames{j},'.PxSE(ii) = std(normPx(ii,:))/sqrt(length(normPx(ii,:)));'])
            eval([phaseNames{j},'.PySE(ii) = std(normPy(ii,:))/sqrt(length(normPy(ii,:)));'])
            eval([phaseNames{j},'.VxSE(ii) = std(normVx(ii,:))/sqrt(length(normVx(ii,:)));'])
            eval([phaseNames{j},'.VySE(ii) = std(normVy(ii,:))/sqrt(length(normVy(ii,:)));'])
            eval([phaseNames{j},'.AxSE(ii) = std(normAx(ii,:))/sqrt(length(normAx(ii,:)));'])
            eval([phaseNames{j},'.AySE(ii) = std(normAy(ii,:))/sqrt(length(normAy(ii,:)));'])
        end
        
        %%% SUBJECT DATA
        eval(['numSubjs = size(',phaseNames{j},'.PxSubj,2);'])
        
        %FOR NORMED SUBJECT MEAN 
        eval(['trLengthSubj = ones(1,numSubjs)*size(',phaseNames{j},'.PxSubjNorm,1);'])
        for jj = 1:numSubjs
            eval(['notLongest = ~isempty(find(isnan(',phaseNames{j},'.PxSubjNorm(:,jj)),1,''first''));'])
            if notLongest    
                eval(['trLengthSubj(jj) = find(isnan(',phaseNames{j},'.PxSubjNorm(:,jj)),1,''first'')-1;'])
            end
        end
%         % Is Floor a bad assumption? Would Rounding up be better or worse?
         avgTrLengthSubj = floor(mean(trLengthSubj));
        
        
        % Resample Subject Trajectories so the same time length
        normSubjPx = zeros(avgTrLengthSubj,numSubjs);
        normSubjPy = zeros(avgTrLengthSubj,numSubjs);
        normSubjVx = zeros(avgTrLengthSubj,numSubjs);
        normSubjVy = zeros(avgTrLengthSubj,numSubjs);
        normSubjAx = zeros(avgTrLengthSubj,numSubjs);
        normSubjAy = zeros(avgTrLengthSubj,numSubjs);
        for jj = 1:numSubjs
            % THIS METHOD USES NORMED MEAN SUBJECT DATA
            eval(['normSubjPx(:,jj) = resample(',phaseNames{j},...
                '.PxSubjNorm(1:trLengthSubj(jj),jj),avgTrLengthSubj,trLengthSubj(jj),0);'])
            eval(['normSubjPy(:,jj) = resample(',phaseNames{j},...
                '.PySubjNorm(1:trLengthSubj(jj),jj),avgTrLengthSubj,trLengthSubj(jj),0);'])            
            eval(['normSubjVx(:,jj) = resample(',phaseNames{j},...
                '.VxSubjNorm(1:trLengthSubj(jj),jj),avgTrLengthSubj,trLengthSubj(jj),0);'])            
            eval(['normSubjVy(:,jj) = resample(',phaseNames{j},...
                '.VySubjNorm(1:trLengthSubj(jj),jj),avgTrLengthSubj,trLengthSubj(jj),0);'])
            eval(['normSubjAx(:,jj) = resample(',phaseNames{j},...
                '.AxSubjNorm(1:trLengthSubj(jj),jj),avgTrLengthSubj,trLengthSubj(jj),0);'])            
            eval(['normSubjAy(:,jj) = resample(',phaseNames{j},...
                '.AySubjNorm(1:trLengthSubj(jj),jj),avgTrLengthSubj,trLengthSubj(jj),0);'])
        end
    
                
        %Save Normalized Trajectories
        eval([phaseNames{j},'.normSubjPx = normSubjPx;'])
        eval([phaseNames{j},'.normSubjPy = normSubjPy;'])
        eval([phaseNames{j},'.normSubjVx = normSubjVx;'])
        eval([phaseNames{j},'.normSubjVy = normSubjVy;'])
        eval([phaseNames{j},'.normSubjAx = normSubjAx;'])
        eval([phaseNames{j},'.normSubjAy = normSubjAy;'])
        
        
        % Find the average trajectory
        for ii = 1:avgTrLengthSubj 
            eval([phaseNames{j},'.PxSubjMean(ii) = mean(normSubjPx(ii,:));'])
            eval([phaseNames{j},'.PySubjMean(ii) = mean(normSubjPy(ii,:));'])
            eval([phaseNames{j},'.VxSubjMean(ii) = mean(normSubjVx(ii,:));'])
            eval([phaseNames{j},'.VySubjMean(ii) = mean(normSubjVy(ii,:));'])
            eval([phaseNames{j},'.AxSubjMean(ii) = mean(normSubjAx(ii,:));'])
            eval([phaseNames{j},'.AySubjMean(ii) = mean(normSubjAy(ii,:));'])            
            eval([phaseNames{j},'.PxSubjSE(ii) = std(normSubjPx(ii,:))/sqrt(length(normSubjPx(ii,:)));'])
            eval([phaseNames{j},'.PySubjSE(ii) = std(normSubjPy(ii,:))/sqrt(length(normSubjPy(ii,:)));'])
            eval([phaseNames{j},'.VxSubjSE(ii) = std(normSubjVx(ii,:))/sqrt(length(normSubjVx(ii,:)));'])
            eval([phaseNames{j},'.VySubjSE(ii) = std(normSubjVy(ii,:))/sqrt(length(normSubjVy(ii,:)));'])
            eval([phaseNames{j},'.AxSubjSE(ii) = std(normSubjAx(ii,:))/sqrt(length(normSubjAx(ii,:)));'])
            eval([phaseNames{j},'.AySubjSE(ii) = std(normSubjAy(ii,:))/sqrt(length(normSubjAy(ii,:)));'])
        end


        % Calculate Standard Error
        eval([phaseNames{j},'.estGainMean = nanmean(',phaseNames{j},'.estGain);'])
        eval([phaseNames{j},'.estGainSE = nanstd(',phaseNames{j},'.estGain)/sqrt(sum(~isnan(',phaseNames{j},'.estGain)));'])

        % Max Perp Mean, STD
        eval([phaseNames{j},'.maxPerpMean = mean(',phaseNames{j},'.maxPerp);'])
        eval([phaseNames{j},'.maxPerpSTD = std(',phaseNames{j},'.maxPerp);'])
        
        % Ang Error Mean, STD
        eval([phaseNames{j},'.angErrMean = mean(',phaseNames{j},'.angErr);'])
        eval([phaseNames{j},'.angErrSTD = std(',phaseNames{j},'.angErr);'])
        else
            eval([phaseNames{j},'.normPx = [];'])
            eval([phaseNames{j},'.normPy = [];'])
            eval([phaseNames{j},'.normVx = [];'])
            eval([phaseNames{j},'.normVy = [];'])
            eval([phaseNames{j},'.normAx = [];'])
            eval([phaseNames{j},'.normAy = [];'])            
            eval([phaseNames{j},'.PxMean = [];'])
            eval([phaseNames{j},'.PyMean = [];'])
            eval([phaseNames{j},'.VxMean = [];'])
            eval([phaseNames{j},'.VyMean = [];'])
            eval([phaseNames{j},'.AxMean = [];'])
            eval([phaseNames{j},'.AyMean = [];'])            
            eval([phaseNames{j},'.PxSE = [];'])
            eval([phaseNames{j},'.PySE = [];'])
            eval([phaseNames{j},'.VxSE = [];'])
            eval([phaseNames{j},'.VySE = [];'])
            eval([phaseNames{j},'.AxSE = [];'])
            eval([phaseNames{j},'.AySE = [];'])
            %
            eval([phaseNames{j},'.normSubjPx = [];'])
            eval([phaseNames{j},'.normSubjPy = [];'])
            eval([phaseNames{j},'.normSubjVx = [];'])
            eval([phaseNames{j},'.normSubjVy = [];'])
            eval([phaseNames{j},'.normSubjAx = [];'])
            eval([phaseNames{j},'.normSubjAy = [];'])            
            eval([phaseNames{j},'.PxSubjMean = [];'])
            eval([phaseNames{j},'.PySubjMean = [];'])
            eval([phaseNames{j},'.VxSubjMean = [];'])
            eval([phaseNames{j},'.VySubjMean = [];'])
            eval([phaseNames{j},'.AxSubjMean = [];'])
            eval([phaseNames{j},'.AySubjMean = [];'])            
            eval([phaseNames{j},'.PxSubjSE = [];'])
            eval([phaseNames{j},'.PySubjSE = [];'])
            eval([phaseNames{j},'.VxSubjSE = [];'])
            eval([phaseNames{j},'.VySubjSE = [];'])
            eval([phaseNames{j},'.AxSubjSE = [];'])
            eval([phaseNames{j},'.AySubjSE = [];'])
            
            eval([phaseNames{j},'.estGainMean = [];'])
            eval([phaseNames{j},'.maxPerpMean = [];'])
            eval([phaseNames{j},'.maxPerpSTD = [];'])
            eval([phaseNames{j},'.angErrMean = [];'])
            eval([phaseNames{j},'.angErrSTD = [];'])
            
        end
    end
    
% FOR NON-TIME NORMALIZED DATA

for ii = 1:length(phaseNames)
    eval(['maxTrialLength = size(',phaseNames{ii},'.Px,1);']);
    eval(['nTrials = size(',phaseNames{ii},'.Px,2);'])
    eval([phaseNames{ii},'.PxMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.PyMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.VxMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.VyMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AxMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AyMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.PxSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.PySENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.VxSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.VySENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AxSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AySENonNorm = zeros(1,maxTrialLength);']);
    for jj = 1:maxTrialLength
        eval([phaseNames{ii},'.PxMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Px(jj,:));'])
        eval([phaseNames{ii},'.PyMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Py(jj,:));'])
        eval([phaseNames{ii},'.VxMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Vx(jj,:));'])
        eval([phaseNames{ii},'.VyMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Vy(jj,:));'])
        eval([phaseNames{ii},'.AxMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Ax(jj,:));'])
        eval([phaseNames{ii},'.AyMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Ay(jj,:));'])
        
        % Commented Out to Inspect STD instead of SE
        eval([phaseNames{ii},'.PxSENonNorm(jj) = nanstd(',phaseNames{ii},'.Px(jj,:))/sqrt(nTrials);'])
        eval([phaseNames{ii},'.PySENonNorm(jj) = nanstd(',phaseNames{ii},'.Py(jj,:))/sqrt(nTrials);'])
        eval([phaseNames{ii},'.VxSENonNorm(jj) = nanstd(',phaseNames{ii},'.Vx(jj,:))/sqrt(nTrials);'])
        eval([phaseNames{ii},'.VySENonNorm(jj) = nanstd(',phaseNames{ii},'.Vy(jj,:))/sqrt(nTrials);'])
        eval([phaseNames{ii},'.AxSENonNorm(jj) = nanstd(',phaseNames{ii},'.Ax(jj,:))/sqrt(nTrials);'])
        eval([phaseNames{ii},'.AySENonNorm(jj) = nanstd(',phaseNames{ii},'.Ay(jj,:))/sqrt(nTrials);'])
    end
    
    % Truncate to mean length.
    eval(['meanLength =  length(',phaseNames{ii},'.PxMean);']);
    eval([phaseNames{ii},'.PxMeanNonNorm = ',phaseNames{ii},'.PxMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.PyMeanNonNorm = ',phaseNames{ii},'.PyMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VxMeanNonNorm = ',phaseNames{ii},'.VxMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VyMeanNonNorm = ',phaseNames{ii},'.VyMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.AxMeanNonNorm = ',phaseNames{ii},'.AxMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.AyMeanNonNorm = ',phaseNames{ii},'.AyMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.PxSENonNorm = ',phaseNames{ii},'.PxSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.PySENonNorm = ',phaseNames{ii},'.PySENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VxSENonNorm = ',phaseNames{ii},'.VxSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VySENonNorm = ',phaseNames{ii},'.VySENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.AxSENonNorm = ',phaseNames{ii},'.AxSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.AySENonNorm = ',phaseNames{ii},'.AySENonNorm(1:meanLength);']);
end

% Different Way to Calculate:
% Average Each Subject's 4 Trials per phase (non-time normalized),
% Then average those averages. SE then becomes the variation between
% averages divided by sqrt(nSubj)
for ii = 1:length(phaseNames)
    eval(['maxTrialLength = size(',phaseNames{ii},'.PxSubj,1);']);
    eval(['nSubj = size(',phaseNames{ii},'.PxSubj,2);'])
    eval([phaseNames{ii},'.PxSubjMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.PySubjMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.VxSubjMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.VySubjMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AxSubjMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AySubjMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.PxSubjSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.PySubjSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.VxSubjSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.VySubjSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AxSubjSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AySubjSENonNorm = zeros(1,maxTrialLength);']);
    for jj = 1:maxTrialLength
        eval([phaseNames{ii},'.PxSubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.PxSubj(jj,:));'])
        eval([phaseNames{ii},'.PySubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.PySubj(jj,:));'])
        eval([phaseNames{ii},'.VxSubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.VxSubj(jj,:));'])
        eval([phaseNames{ii},'.VySubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.VySubj(jj,:));'])
        eval([phaseNames{ii},'.AxSubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.AxSubj(jj,:));'])
        eval([phaseNames{ii},'.AySubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.AySubj(jj,:));'])
        
        % Commented Out to Inspect STD instead of SE
        eval([phaseNames{ii},'.PxSubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.PxSubj(jj,:))/sqrt(nSubj);'])
        eval([phaseNames{ii},'.PySubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.PySubj(jj,:))/sqrt(nSubj);'])
        eval([phaseNames{ii},'.VxSubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.VxSubj(jj,:))/sqrt(nSubj);'])
        eval([phaseNames{ii},'.VySubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.VySubj(jj,:))/sqrt(nSubj);'])
        eval([phaseNames{ii},'.AxSubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.AxSubj(jj,:))/sqrt(nSubj);'])
        eval([phaseNames{ii},'.AySubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.AySubj(jj,:))/sqrt(nSubj);'])
    end
        
    if age == 'Y' || age == 'y'
        trialLengths = [186,98,140,131,98,103,135,132];
    else
        trialLengths = [152,127,157,127,140,115,158,117]; %TBD
    end
    meanLength = trialLengths(ii);

    
    eval([phaseNames{ii},'.PxSubjMeanNonNorm = ',phaseNames{ii},'.PxSubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.PySubjMeanNonNorm = ',phaseNames{ii},'.PySubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VxSubjMeanNonNorm = ',phaseNames{ii},'.VxSubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VySubjMeanNonNorm = ',phaseNames{ii},'.VySubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.AxSubjMeanNonNorm = ',phaseNames{ii},'.AxSubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.AySubjMeanNonNorm = ',phaseNames{ii},'.AySubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.PxSubjSENonNorm = ',phaseNames{ii},'.PxSubjSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.PySubjSENonNorm = ',phaseNames{ii},'.PySubjSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VxSubjSENonNorm = ',phaseNames{ii},'.VxSubjSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VySubjSENonNorm = ',phaseNames{ii},'.VySubjSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.AxSubjSENonNorm = ',phaseNames{ii},'.AxSubjSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.AySubjSENonNorm = ',phaseNames{ii},'.AySubjSENonNorm(1:meanLength);']);
    
end

% Add Curl Gain Information
for j = 1:length(phaseNames)
    if j>=3 && j <=6
        eval([phaseNames{j},'.curlGain = -20;'])
    else
        eval([phaseNames{j},'.curlGain = 0;'])
    end
end
    
% Combine into Structure for Variable Out
for j = 1:length(phaseNames)
    eval(['D.',phaseNames{j},' = ',phaseNames{j},';'])  
end