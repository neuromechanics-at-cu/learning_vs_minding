%% Process Channel Trial Data from Old or Young Subects
% 
% Throws out bad trials/unnecessary data
% Processes averages for each subject group and phase of the experiment

function D = processRawChannelTrialData(age,trial_bin_type)
% 'subj' is 3-letter acronym, ex: 'FTN'
% 'age' is age group, ex: 'Y' or 'O'

addpath('../../../Scripts/Supporting Functions')

% curDir = pwd;
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
    
phaseNames = {'EN1','LN1','EF1','LF1','EF2','LF2','EN2','LN2'};

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

    
    % Get Trial Numbers for Specified Phase
    for j = 1:length(phaseNames)

        eval(['ctrialNums.',phaseNames{j},' = intersect(allCTrialNums,phaseTrialNums.',phaseNames{j},');'])
        
        % Initialize Velocity Filter Trial Structure
        eval([phaseNames{j},'.ctrialNums{k}.tooLong = [];'])
        eval([phaseNames{j},'.ctrialNums{k}.tooShort = [];'])
        eval([phaseNames{j},'.ctrialNums{k}.justRight = [];'])

        % Find Number of Trials
        eval(['nTrials = length(ctrialNums.',phaseNames{j},');'])

        % Find Max Trial Length
        maxLength = 0;
        for m = 1:nTrials
            eval(['lengths(m) = length(T.framedata(1,ctrialNums.',phaseNames{j},'(m)).x);'])
            if lengths(m) > maxLength
               maxLength = lengths(m); 
            end
        end
        
        % Add additional trials to TR structure
        % Preallocate Matrix
        eval(['TR.ctrials.Px.',phaseNames{j},'.out = nan(maxLength,nTrials);'])
        eval(['TR.ctrials.Py.',phaseNames{j},'.out = nan(maxLength,nTrials);'])
        eval(['TR.ctrials.Vx.',phaseNames{j},'.out = nan(maxLength,nTrials);'])
        eval(['TR.ctrials.Vy.',phaseNames{j},'.out = nan(maxLength,nTrials);'])

        % Fill Matrix with Trial Data
        for m = 1:nTrials
            eval(['TR.ctrials.Px.',phaseNames{j},'.out(1:lengths(m),m) = T.framedata(1,ctrialNums.',phaseNames{j},'(m)).x;'])
            eval(['TR.ctrials.Py.',phaseNames{j},'.out(1:lengths(m),m) = T.framedata(1,ctrialNums.',phaseNames{j},'(m)).y;'])
            eval(['TR.ctrials.Vx.',phaseNames{j},'.out(1:lengths(m),m) = T.framedata(1,ctrialNums.',phaseNames{j},'(m)).vx;'])
            eval(['TR.ctrials.Vy.',phaseNames{j},'.out(1:lengths(m),m) = T.framedata(1,ctrialNums.',phaseNames{j},'(m)).vy;'])
            eval(['TR.ctrials.Fx.',phaseNames{j},'.out(1:lengths(m),m) = T.framedata(1,ctrialNums.',phaseNames{j},'(m)).fx;'])
            eval(['TR.ctrials.Fy.',phaseNames{j},'.out(1:lengths(m),m) = T.framedata(1,ctrialNums.',phaseNames{j},'(m)).fy;'])
        end
        
        % Add Ax, Ay
        eval(['TR.ctrials.Ax.',phaseNames{j},'.out = nan(maxLength,nTrials);'])
        eval(['TR.ctrials.Ay.',phaseNames{j},'.out = nan(maxLength,nTrials);'])
        for m = 1:nTrials
            eval(['TR.ctrials.Ax.',phaseNames{j},'.out(1:lengths(m),m) = NumDiff_3pt(T.framedata(1,ctrialNums.',phaseNames{j},'(m)).vx,[1/200:1/200:lengths(m)/200]'');'])
            eval(['TR.ctrials.Ay.',phaseNames{j},'.out(1:lengths(m),m) = NumDiff_3pt(T.framedata(1,ctrialNums.',phaseNames{j},'(m)).vy,[1/200:1/200:lengths(m)/200]'');'])
            % Add butterworth, zero shift filter, to smooth acceleration
            [bb,aa] = butter(9,.1,'low'); % 0.2 ~ 20 Hz
            eval(['TR.ctrials.Ax.',phaseNames{j},'.out(1:lengths(m),m) = filtfilt(bb,aa,TR.ctrials.Ax.',phaseNames{j},'.out(1:lengths(m),m));'])
            eval(['TR.ctrials.Ay.',phaseNames{j},'.out(1:lengths(m),m) = filtfilt(bb,aa,TR.ctrials.Ay.',phaseNames{j},'.out(1:lengths(m),m));'])
        end
            
            
        
        if nTrials >= 1
        maxLength = 0;
        % Find Movement Start (Vy > 0.03 m/s) and Movement End (Vy < 0.03 m/s)   
        for m = 1:nTrials 
%             eval(['start = find(TR.ctrials.Vy.',phaseNames{j},...
%                 '.out(:,m) > 0.03 & (abs(TR.ctrials.Px.',phaseNames{j},...
%                 '.out(:,m)) > 0.008 | abs(TR.ctrials.Py.',phaseNames{j},...
%                 '.out(:,m) + 0.1) > 0.008),1);'])
%             if start>15
%                 start = start - 15;
%             else 
%                 start = 1;
%             end
            eval(['start = find(TR.ctrials.Vy.',phaseNames{j},...
                '.out(:,m) > 0.03,1);'])

            eval(['finish = start - 2 + ',...
                'find(TR.ctrials.Vy.',phaseNames{j},'.out(start:end,m) < 0.03',...
                '& abs(TR.ctrials.Px.',phaseNames{j},'.out(start:end,m)) < 0.008',...
                '& abs(TR.ctrials.Py.',phaseNames{j},'.out(start:end,m) - 0.1) < 0.008,1);'])
            % CMH Added 5/15/22 after invesitgating Fx/Vy plots in LF2
            % Loosen target bounds -- anything on second half of the reach
            eval(['vel_finish = start - 2 + find(TR.ctrials.Vy.',phaseNames{j},'.out(start:end,m) < 0.03',...
                                                '& TR.ctrials.Py.',phaseNames{j},'.out(start:end,m) > 0,1);'])
            if isempty(finish) || ((finish - vel_finish) >= 50)
                finish = vel_finish;
% %                 finish
%                 if j == 2
%                     finish
%                     figure(1)
%                     hold on
%                     plot(TR.ctrials.Vy.LN1.out(start:finish))
%                 end
            end
            if isempty(finish)
                eval(['finish = length(TR.ctrials.Vy.',phaseNames{j},'.out(:,m));'])
            end
            % Resave Trajectories to New Length
            eval(['TR.ctrials.Px.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.ctrials.Px.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.ctrials.Py.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.ctrials.Py.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.ctrials.Vx.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.ctrials.Vx.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.ctrials.Vy.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.ctrials.Vy.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.ctrials.Fx.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.ctrials.Fx.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.ctrials.Fy.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.ctrials.Fy.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.ctrials.Ax.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.ctrials.Ax.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.ctrials.Ay.',phaseNames{j},'.out(1:(finish-start+1),m) = TR.ctrials.Ay.',phaseNames{j},'.out(start:finish,m);'])
            eval(['TR.ctrials.Px.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
            eval(['TR.ctrials.Py.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
            eval(['TR.ctrials.Vx.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
            eval(['TR.ctrials.Vy.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
            eval(['TR.ctrials.Fx.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
            eval(['TR.ctrials.Fy.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
            eval(['TR.ctrials.Ax.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
            eval(['TR.ctrials.Ay.',phaseNames{j},'.out(finish-start+2:end,m) = NaN;'])
        end
        % Remove NaNs
        eval(['[inan,jnan] = find(~isnan(TR.ctrials.Px.',phaseNames{j},'.out));'])
        maxLength = max(inan);
        eval(['TR.ctrials.Px.',phaseNames{j},'.out = TR.ctrials.Px.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Py.',phaseNames{j},'.out = TR.ctrials.Py.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Vx.',phaseNames{j},'.out = TR.ctrials.Vx.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Vy.',phaseNames{j},'.out = TR.ctrials.Vy.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Fx.',phaseNames{j},'.out = TR.ctrials.Fx.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Fy.',phaseNames{j},'.out = TR.ctrials.Fy.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Ax.',phaseNames{j},'.out = TR.ctrials.Ax.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Ay.',phaseNames{j},'.out = TR.ctrials.Ay.',phaseNames{j},'.out(1:maxLength,:);'])

        end
        
        % Eliminate Botched Trials 
        % (Ex: Subj-'FTN', Phase-EF2)
         minLength = 50;%50     60  (300 ms)
         maxLength = 300;%300   120 (600 ms)
         if eval(['size(TR.ctrials.Px.',phaseNames{j},'.out,2)==0'])
                 eval(['ctrialNums.',phaseNames{j},' = [];']) 
                 eval(['TR.ctrials.Px.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Py.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Vx.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Vy.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Fx.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Fy.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Ax.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Ay.',phaseNames{j},'.out = [];'])             
         elseif eval(['size(TR.ctrials.Px.',phaseNames{j},'.out,2)==1'])
             if eval(['length(TR.ctrials.Px.',phaseNames{j},'.out) < minLength'])
                 % CMHNEW
                 eval([phaseNames{j},'.ctrialNums{k}.tooShort = [',phaseNames{j},'.ctrialNums{k}.tooShort,ctrialNums.',phaseNames{j},'];'])
                 %
                 eval(['ctrialNums.',phaseNames{j},' = [];']) 
                 eval(['TR.ctrials.Px.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Py.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Vx.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Vy.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Fx.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Fy.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Ax.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Ay.',phaseNames{j},'.out = [];'])
             elseif eval(['length(TR.ctrials.Px.',phaseNames{j},'.out) > maxLength'])
                 % CMHNEW
                 eval([phaseNames{j},'.ctrialNums{k}.tooLong = [',phaseNames{j},'.ctrialNums{k}.tooLong,ctrialNums.',phaseNames{j},'];'])
                 %
                 eval(['ctrialNums.',phaseNames{j},' = [];']) 
                 eval(['TR.ctrials.Px.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Py.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Vx.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Vy.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Fx.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Fy.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Ax.',phaseNames{j},'.out = [];'])
                 eval(['TR.ctrials.Ay.',phaseNames{j},'.out = [];'])
             end
             % CMHNEW
             eval([phaseNames{j},'.ctrialNums{k}.justRight = [',phaseNames{j},'.ctrialNums{k}.justRight,ctrialNums.',phaseNames{j},'];'])
             %
             
         else
             % Remove if Too Short
             eval(['[II,JJ] = find(isnan(TR.ctrials.Px.',phaseNames{j},'.out));'])
             III = find(II < minLength);
             % CMHNEW
             eval([phaseNames{j},'.ctrialNums{k}.tooShort = [',phaseNames{j},'.ctrialNums{k}.tooShort,unique(JJ(III))];'])
             %
             eval(['ctrialNums.',phaseNames{j},'(JJ(III)) = [];']) 
             eval(['TR.ctrials.Px.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             eval(['TR.ctrials.Py.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             eval(['TR.ctrials.Vx.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             eval(['TR.ctrials.Vy.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             eval(['TR.ctrials.Fx.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             eval(['TR.ctrials.Fy.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             eval(['TR.ctrials.Ax.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             eval(['TR.ctrials.Ay.',phaseNames{j},'.out(:,JJ(III)) = [];'])
             % Remove if Too Long
             eval(['[nII,nJJ] = find(~isnan(TR.ctrials.Px.',phaseNames{j},'.out));'])
             nIII = find(nII > maxLength);
             % CMHNEW
             eval([phaseNames{j},'.ctrialNums{k}.tooLong = [',phaseNames{j},'.ctrialNums{k}.tooLong,unique(nJJ(nIII))];'])
             %
             eval(['ctrialNums.',phaseNames{j},'(nJJ(nIII)) = [];'])
             eval(['TR.ctrials.Px.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             eval(['TR.ctrials.Py.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             eval(['TR.ctrials.Vx.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             eval(['TR.ctrials.Vy.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             eval(['TR.ctrials.Fx.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             eval(['TR.ctrials.Fy.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             eval(['TR.ctrials.Ax.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             eval(['TR.ctrials.Ay.',phaseNames{j},'.out(:,nJJ(nIII)) = [];'])
             % CMHNEW
             eval([phaseNames{j},'.ctrialNums{k}.justRight = [',phaseNames{j},'.ctrialNums{k}.justRight,ctrialNums.',phaseNames{j},'];'])
             %
         end
         
        % Re-remove NaNs
        eval(['[inan,jnan] = find(~isnan(TR.ctrials.Px.',phaseNames{j},'.out));'])
        maxLength = max(inan);
        eval(['TR.ctrials.Px.',phaseNames{j},'.out = TR.ctrials.Px.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Py.',phaseNames{j},'.out = TR.ctrials.Py.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Vx.',phaseNames{j},'.out = TR.ctrials.Vx.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Vy.',phaseNames{j},'.out = TR.ctrials.Vy.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Fx.',phaseNames{j},'.out = TR.ctrials.Fx.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Fy.',phaseNames{j},'.out = TR.ctrials.Fy.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Ax.',phaseNames{j},'.out = TR.ctrials.Ax.',phaseNames{j},'.out(1:maxLength,:);'])
        eval(['TR.ctrials.Ay.',phaseNames{j},'.out = TR.ctrials.Ay.',phaseNames{j},'.out(1:maxLength,:);'])
        
        if isempty(maxLength)
                eval(['TR.ctrials.Px.',phaseNames{j},'.mean = [];'])
                eval(['TR.ctrials.Py.',phaseNames{j},'.mean = [];'])
                eval(['TR.ctrials.Vx.',phaseNames{j},'.mean = [];'])
                eval(['TR.ctrials.Vy.',phaseNames{j},'.mean = [];'])
                eval(['TR.ctrials.Fx.',phaseNames{j},'.mean = [];'])
                eval(['TR.ctrials.Fy.',phaseNames{j},'.mean = [];'])
                eval(['TR.ctrials.Ax.',phaseNames{j},'.mean = [];'])
                eval(['TR.ctrials.Ay.',phaseNames{j},'.mean = [];'])
                eval(['TR.ctrials.Px.',phaseNames{j},'.std = [];'])
                eval(['TR.ctrials.Py.',phaseNames{j},'.std = [];'])
                eval(['TR.ctrials.Vx.',phaseNames{j},'.std = [];'])
                eval(['TR.ctrials.Vy.',phaseNames{j},'.std = [];'])
                eval(['TR.ctrials.Fx.',phaseNames{j},'.std = [];'])
                eval(['TR.ctrials.Fy.',phaseNames{j},'.std = [];'])
                eval(['TR.ctrials.Ax.',phaseNames{j},'.std = [];'])
                eval(['TR.ctrials.Ay.',phaseNames{j},'.std = [];'])
        else
            % Find Within Subject Mean, STD
            for p = 1:maxLength 
                eval(['TR.ctrials.Px.',phaseNames{j},'.mean(p,1) = nanmean(TR.ctrials.Px.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Py.',phaseNames{j},'.mean(p,1) = nanmean(TR.ctrials.Py.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Vx.',phaseNames{j},'.mean(p,1) = nanmean(TR.ctrials.Vx.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Vy.',phaseNames{j},'.mean(p,1) = nanmean(TR.ctrials.Vy.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Fx.',phaseNames{j},'.mean(p,1) = nanmean(TR.ctrials.Fx.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Fy.',phaseNames{j},'.mean(p,1) = nanmean(TR.ctrials.Fy.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Ax.',phaseNames{j},'.mean(p,1) = nanmean(TR.ctrials.Ax.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Ay.',phaseNames{j},'.mean(p,1) = nanmean(TR.ctrials.Ay.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Px.',phaseNames{j},'.std(p,1) = nanstd(TR.ctrials.Px.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Py.',phaseNames{j},'.std(p,1) = nanstd(TR.ctrials.Py.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Vx.',phaseNames{j},'.std(p,1) = nanstd(TR.ctrials.Vx.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Vy.',phaseNames{j},'.std(p,1) = nanstd(TR.ctrials.Vy.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Fx.',phaseNames{j},'.std(p,1) = nanstd(TR.ctrials.Fx.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Fy.',phaseNames{j},'.std(p,1) = nanstd(TR.ctrials.Fy.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Ax.',phaseNames{j},'.std(p,1) = nanstd(TR.ctrials.Ax.',phaseNames{j},'.out(p,:));'])
                eval(['TR.ctrials.Ay.',phaseNames{j},'.std(p,1) = nanstd(TR.ctrials.Ay.',phaseNames{j},'.out(p,:));'])
            end 
        end
        
        % Re-Find Number of Trials
        eval(['nTrials = length(ctrialNums.',phaseNames{j},');'])
        
        % Track Max Subject Trial Length
        eval([phaseNames{j},'.SubjTrLength(k) = length(TR.ctrials.Px.',phaseNames{j},'.mean);'])
        % Track Mean Subject Trial Length
        eval(['trLength = ones(1,nTrials)*size(TR.ctrials.Px.',phaseNames{j},'.out,1);'])
        for jj = 1:nTrials
            eval(['notLongest = ~isempty(find(isnan(TR.ctrials.Px.',phaseNames{j},'.out(:,jj)),1,''first''));'])
            if notLongest    
                eval(['trLength(jj) = find(isnan(TR.ctrials.Px.',phaseNames{j},'.out(:,jj)),1,''first'')-1;'])
            end
        end
        % Not going to use a floor call yet.
        eval([phaseNames{j},'.SubjMeanTrLength(k) = mean(trLength);'])

        if nTrials == 0
            avgTrLengthSubj = floor(mean(trLength));
            PxSubjNorm = [];
            PySubjNorm = [];
            VxSubjNorm = [];
            VySubjNorm = [];
            FxSubjNorm = [];
            FySubjNorm = [];
            AxSubjNorm = [];
            AySubjNorm = [];
            eval(['TR.ctrials.Px.',phaseNames{j},'.normed = [];'])
            eval(['TR.ctrials.Py.',phaseNames{j},'.normed = [];'])
            eval(['TR.ctrials.Vx.',phaseNames{j},'.normed = [];'])
            eval(['TR.ctrials.Vy.',phaseNames{j},'.normed = [];'])
            eval(['TR.ctrials.Fx.',phaseNames{j},'.normed = [];'])
            eval(['TR.ctrials.Fy.',phaseNames{j},'.normed = [];'])
            eval(['TR.ctrials.Ax.',phaseNames{j},'.normed = [];'])
            eval(['TR.ctrials.Ay.',phaseNames{j},'.normed = [];'])
        else
            % Re-sample individual trials to subject mean length, then average
            % for subject mean-normed trajectory
            % Resample Subject Trajectories so the same time length
            avgTrLengthSubj = floor(mean(trLength));
            PxSubjNorm = zeros(avgTrLengthSubj,nTrials);
            PySubjNorm = zeros(avgTrLengthSubj,nTrials);
            VxSubjNorm = zeros(avgTrLengthSubj,nTrials);
            VySubjNorm = zeros(avgTrLengthSubj,nTrials);
            FxSubjNorm = zeros(avgTrLengthSubj,nTrials);
            FySubjNorm = zeros(avgTrLengthSubj,nTrials);
            AxSubjNorm = zeros(avgTrLengthSubj,nTrials);
            AySubjNorm = zeros(avgTrLengthSubj,nTrials);
            for jj = 1:nTrials
                eval(['PxSubjNorm(:,jj) = resample(TR.ctrials.Px.',phaseNames{j},...
                    '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])
                eval(['PySubjNorm(:,jj) = resample(TR.ctrials.Py.',phaseNames{j},...
                    '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])            
                eval(['VxSubjNorm(:,jj) = resample(TR.ctrials.Vx.',phaseNames{j},...
                    '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])            
                eval(['VySubjNorm(:,jj) = resample(TR.ctrials.Vy.',phaseNames{j},...
                    '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])
                eval(['FxSubjNorm(:,jj) = resample(TR.ctrials.Fx.',phaseNames{j},...
                    '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])            
                eval(['FySubjNorm(:,jj) = resample(TR.ctrials.Fy.',phaseNames{j},...
                    '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])
                eval(['AxSubjNorm(:,jj) = resample(TR.ctrials.Ax.',phaseNames{j},...
                    '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])            
                eval(['AySubjNorm(:,jj) = resample(TR.ctrials.Ay.',phaseNames{j},...
                    '.out(1:trLength(jj),jj),avgTrLengthSubj,trLength(jj),0);'])
            end
            for ii = 1:avgTrLengthSubj
                eval(['TR.ctrials.Px.',phaseNames{j},'.normed(ii,:) = mean(PxSubjNorm(ii,:));'])
                eval(['TR.ctrials.Py.',phaseNames{j},'.normed(ii,:) = mean(PySubjNorm(ii,:));'])
                eval(['TR.ctrials.Vx.',phaseNames{j},'.normed(ii,:) = mean(VxSubjNorm(ii,:));'])
                eval(['TR.ctrials.Vy.',phaseNames{j},'.normed(ii,:) = mean(VySubjNorm(ii,:));'])
                eval(['TR.ctrials.Fx.',phaseNames{j},'.normed(ii,:) = mean(FxSubjNorm(ii,:));'])
                eval(['TR.ctrials.Fy.',phaseNames{j},'.normed(ii,:) = mean(FySubjNorm(ii,:));'])
                eval(['TR.ctrials.Ax.',phaseNames{j},'.normed(ii,:) = mean(AxSubjNorm(ii,:));'])
                eval(['TR.ctrials.Ay.',phaseNames{j},'.normed(ii,:) = mean(AySubjNorm(ii,:));'])
            end
        end
        
    % Assemble Subject Mean/STD into Structure
        if k == 1
            eval([phaseNames{j},'.PxSubj = TR.ctrials.Px.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.PySubj = TR.ctrials.Py.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.VxSubj = TR.ctrials.Vx.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.VySubj = TR.ctrials.Vy.',phaseNames{j},'.mean;']) 
            eval([phaseNames{j},'.FxSubj = TR.ctrials.Fx.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.FySubj = TR.ctrials.Fy.',phaseNames{j},'.mean;']) 
            eval([phaseNames{j},'.AxSubj = TR.ctrials.Ax.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.AySubj = TR.ctrials.Ay.',phaseNames{j},'.mean;']) 
        elseif nTrials > 0
            eval(['trialLength = size(TR.ctrials.Px.',phaseNames{j},'.mean,1);'])
            eval(['matrixLength = size(',phaseNames{j},'.PxSubj,1);'])
            if trialLength > matrixLength
                eval([phaseNames{j},'.PxSubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.PxSubj,2));'])
                eval([phaseNames{j},'.PySubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.PySubj,2));'])
                eval([phaseNames{j},'.VxSubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.VxSubj,2));'])
                eval([phaseNames{j},'.VySubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.VySubj,2));'])  
                eval([phaseNames{j},'.FxSubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.FxSubj,2));'])
                eval([phaseNames{j},'.FySubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.FySubj,2));'])    
                eval([phaseNames{j},'.AxSubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.AxSubj,2));'])
                eval([phaseNames{j},'.AySubj(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.AySubj,2));'])                  
            elseif trialLength < matrixLength
                eval(['TR.ctrials.Px.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Py.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Vx.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Vy.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Fx.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Fy.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Ax.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Ay.',phaseNames{j},'.mean(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
            end
            eval([phaseNames{j},'.PxSubj(:,end+1) = TR.ctrials.Px.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.PySubj(:,end+1) = TR.ctrials.Py.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.VxSubj(:,end+1) = TR.ctrials.Vx.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.VySubj(:,end+1) = TR.ctrials.Vy.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.FxSubj(:,end+1) = TR.ctrials.Fx.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.FySubj(:,end+1) = TR.ctrials.Fy.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.AxSubj(:,end+1) = TR.ctrials.Ax.',phaseNames{j},'.mean;'])
            eval([phaseNames{j},'.AySubj(:,end+1) = TR.ctrials.Ay.',phaseNames{j},'.mean;'])
        end

        % Assemble Subject NORMED Mean/STD into Structure
        if k == 1 
            eval([phaseNames{j},'.PxSubjNorm = TR.ctrials.Px.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.PySubjNorm = TR.ctrials.Py.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.VxSubjNorm = TR.ctrials.Vx.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.VySubjNorm = TR.ctrials.Vy.',phaseNames{j},'.normed;']) 
            eval([phaseNames{j},'.FxSubjNorm = TR.ctrials.Fx.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.FySubjNorm = TR.ctrials.Fy.',phaseNames{j},'.normed;']) 
            eval([phaseNames{j},'.AxSubjNorm = TR.ctrials.Ax.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.AySubjNorm = TR.ctrials.Ay.',phaseNames{j},'.normed;']) 
        elseif nTrials > 0
            eval(['trialLength = size(TR.ctrials.Px.',phaseNames{j},'.normed,1);'])
            eval(['matrixLength = size(',phaseNames{j},'.PxSubjNorm,1);'])
            if trialLength > matrixLength
                eval([phaseNames{j},'.PxSubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.PxSubjNorm,2));'])
                eval([phaseNames{j},'.PySubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.PySubjNorm,2));'])
                eval([phaseNames{j},'.VxSubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.VxSubjNorm,2));'])
                eval([phaseNames{j},'.VySubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.VySubjNorm,2));'])  
                eval([phaseNames{j},'.FxSubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.FxSubjNorm,2));'])
                eval([phaseNames{j},'.FySubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.FySubjNorm,2));'])  
                eval([phaseNames{j},'.AxSubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.AxSubjNorm,2));'])
                eval([phaseNames{j},'.AySubjNorm(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.AySubjNorm,2));'])                  
            elseif trialLength < matrixLength
                eval(['TR.ctrials.Px.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Py.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Vx.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Vy.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Fx.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Fy.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Ax.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
                eval(['TR.ctrials.Ay.',phaseNames{j},'.normed(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,1);'])
            end
            eval([phaseNames{j},'.PxSubjNorm(:,end+1) = TR.ctrials.Px.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.PySubjNorm(:,end+1) = TR.ctrials.Py.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.VxSubjNorm(:,end+1) = TR.ctrials.Vx.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.VySubjNorm(:,end+1) = TR.ctrials.Vy.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.FxSubjNorm(:,end+1) = TR.ctrials.Fx.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.FySubjNorm(:,end+1) = TR.ctrials.Fy.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.AxSubjNorm(:,end+1) = TR.ctrials.Ax.',phaseNames{j},'.normed;'])
            eval([phaseNames{j},'.AySubjNorm(:,end+1) = TR.ctrials.Ay.',phaseNames{j},'.normed;'])
        end

        if k == 1
            eval([phaseNames{j},'.Px = TR.ctrials.Px.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Py = TR.ctrials.Py.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Vx = TR.ctrials.Vx.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Vy = TR.ctrials.Vy.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Fx = TR.ctrials.Fx.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Fy = TR.ctrials.Fy.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Ax = TR.ctrials.Ax.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Ay = TR.ctrials.Ay.',phaseNames{j},'.out;'])
        elseif nTrials > 0
            eval(['trialLength = size(TR.ctrials.Px.',phaseNames{j},'.out,1);'])
            eval(['matrixLength = size(',phaseNames{j},'.Px,1);'])
            if trialLength > matrixLength
                eval([phaseNames{j},'.Px(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Px,2));'])
                eval([phaseNames{j},'.Py(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Py,2));'])
                eval([phaseNames{j},'.Vx(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Vx,2));'])
                eval([phaseNames{j},'.Vy(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Vy,2));'])            
                eval([phaseNames{j},'.Fx(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Fx,2));'])
                eval([phaseNames{j},'.Fy(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Fy,2));']) 
                eval([phaseNames{j},'.Ax(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Ax,2));'])
                eval([phaseNames{j},'.Ay(matrixLength+1:trialLength,:) = nan(trialLength - matrixLength,size(',phaseNames{j},'.Ay,2));'])            
            elseif trialLength < matrixLength
                eval(['TR.ctrials.Px.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
                eval(['TR.ctrials.Py.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
                eval(['TR.ctrials.Vx.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
                eval(['TR.ctrials.Vy.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
                eval(['TR.ctrials.Fx.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
                eval(['TR.ctrials.Fy.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
                eval(['TR.ctrials.Ax.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
                eval(['TR.ctrials.Ay.',phaseNames{j},'.out(trialLength+1:matrixLength,:) = nan(matrixLength - trialLength,nTrials);'])
            end
            eval([phaseNames{j},'.Px(:,end+1:end+nTrials) = TR.ctrials.Px.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Py(:,end+1:end+nTrials) = TR.ctrials.Py.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Vx(:,end+1:end+nTrials) = TR.ctrials.Vx.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Vy(:,end+1:end+nTrials) = TR.ctrials.Vy.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Fx(:,end+1:end+nTrials) = TR.ctrials.Fx.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Fy(:,end+1:end+nTrials) = TR.ctrials.Fy.',phaseNames{j},'.out;'])    
            eval([phaseNames{j},'.Ax(:,end+1:end+nTrials) = TR.ctrials.Ax.',phaseNames{j},'.out;'])
            eval([phaseNames{j},'.Ay(:,end+1:end+nTrials) = TR.ctrials.Ay.',phaseNames{j},'.out;'])            
        end
        
        % Get Curl Field Estimation

        %      Calculates average if phase has multiple catch trials
        estGain = 0;
        FxMax = 0;
        eval(['nTemp = length(ctrialNums.',phaseNames{j},');']);
        for i = 1:nTemp
            eval(['estGain = estGain + fitgain(T.framedata(1,ctrialNums.',phaseNames{j},'(i)).fx,T.framedata(1,ctrialNums.',phaseNames{j},'(i)).vy);'])
            eval(['[~,t_ind] = max(abs(T.framedata(1,ctrialNums.',phaseNames{j},'(i)).fx));'])
            eval(['FxMax = FxMax + T.framedata(1,ctrialNums.',phaseNames{j},'(i)).fx(t_ind);'])
        end
%         eval([phaseNames{j},'.estGain(k) = -1*estGain/length(ctrialNums.',phaseNames{j},');'])%Takes averages, converts to negative (from abs. value)
        eval([phaseNames{j},'.estGain(k) = estGain/length(ctrialNums.',phaseNames{j},');'])%Takes averages,
        eval([phaseNames{j},'.FxMax(k) = FxMax/length(ctrialNums.',phaseNames{j},');'])%Takes averages,

        %Max Perp and Angular Error
        if k == 1
           eval([phaseNames{j},'.maxPerp = [];'])
            if nTrials > 0
                for i = 1:nTrials
                    eval(['[~,maxPerpInd] = max(abs(TR.ctrials.Px.',phaseNames{j},'.out(:,i)));'])
                    eval([phaseNames{j},'.maxPerp = TR.ctrials.Px.',phaseNames{j},'.out(maxPerpInd,i);'])
                    
                    eval(['[~,maxVelInd] = max(TR.ctrials.Vx.',phaseNames{j},...
                        '.out(:,i).^2 + TR.ctrials.Vy.',phaseNames{j},'.out(:,i).^2);'])
                    eval([phaseNames{j},'.angErr = acosd(dot([TR.ctrials.Px.',phaseNames{j},...
                        '.out(maxVelInd,i) - 0; TR.ctrials.Py.',phaseNames{j},'.out(maxVelInd,i) + 0.1],[0;1])/(norm([TR.ctrials.Px.',phaseNames{j},...
                        '.out(maxVelInd,i) - 0; TR.ctrials.Py.',phaseNames{j},'.out(maxVelInd,i) + 0.1])*norm([0;1])));'])
                end
            end
        else
            eval(['totTrials = length(',phaseNames{j},'.maxPerp);'])
            for i = 1:nTrials
                eval(['[~,maxPerpInd] = max(abs(TR.ctrials.Px.',phaseNames{j},'.out(:,i)));'])
                eval([phaseNames{j},'.maxPerp(totTrials+i) = TR.ctrials.Px.',phaseNames{j},'.out(maxPerpInd,i);'])
                
                eval(['[~,maxVelInd] = max(TR.ctrials.Vx.',phaseNames{j},...
                        '.out(:,i).^2 + TR.ctrials.Vy.',phaseNames{j},'.out(:,i).^2);'])
                eval([phaseNames{j},'.angErr(totTrials+i) = acosd(dot([TR.ctrials.Px.',phaseNames{j},...
                        '.out(maxVelInd,i) - 0; TR.ctrials.Py.',phaseNames{j},'.out(maxVelInd,i) + 0.1],[0;1])/(norm([TR.ctrials.Px.',phaseNames{j},...
                        '.out(maxVelInd,i) - 0; TR.ctrials.Py.',phaseNames{j},'.out(maxVelInd,i) + 0.1])*norm([0;1])));'])

            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Average Subject Data

% if length(subjs) > 1
    for j = 1:length(phaseNames)
        % Find Valid Trials for Given Phase
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
            normFx = zeros(avgTrLength,numTrials);
            normFy = zeros(avgTrLength,numTrials);
            normAx = zeros(avgTrLength,numTrials);
            normAy = zeros(avgTrLength,numTrials);
            for jj = 1:numTrials
                eval(['normPx(:,jj) = resample(',phaseNames{j},...
                    '.Px(1:trLength(jj),jj),avgTrLength,trLength(jj),0);'])
                eval(['normPy(:,jj) = resample(',phaseNames{j},...
                    '.Py(1:trLength(jj),jj),avgTrLength,trLength(jj),0);'])            
                eval(['normVx(:,jj) = resample(',phaseNames{j},...
                    '.Vx(1:trLength(jj),jj),avgTrLength,trLength(jj),0);'])            
                eval(['normVy(:,jj) = resample(',phaseNames{j},...
                    '.Vy(1:trLength(jj),jj),avgTrLength,trLength(jj),0);'])
                eval(['normFx(:,jj) = resample(',phaseNames{j},...
                    '.Fx(1:trLength(jj),jj),avgTrLength,trLength(jj),0);'])            
                eval(['normFy(:,jj) = resample(',phaseNames{j},...
                    '.Fy(1:trLength(jj),jj),avgTrLength,trLength(jj),0);'])
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
            eval([phaseNames{j},'.normFx = normFx;'])
            eval([phaseNames{j},'.normFy = normFy;'])
            eval([phaseNames{j},'.normAx = normAx;'])
            eval([phaseNames{j},'.normAy = normAy;'])
        
            % Find the average trajectory
            for ii = 1:avgTrLength
                eval([phaseNames{j},'.PxMean(ii) = mean(normPx(ii,:));'])
                eval([phaseNames{j},'.PyMean(ii) = mean(normPy(ii,:));'])
                eval([phaseNames{j},'.VxMean(ii) = mean(normVx(ii,:));'])
                eval([phaseNames{j},'.VyMean(ii) = mean(normVy(ii,:));'])
                eval([phaseNames{j},'.FxMean(ii) = mean(normFx(ii,:));'])
                eval([phaseNames{j},'.FyMean(ii) = mean(normFy(ii,:));'])       
                eval([phaseNames{j},'.AxMean(ii) = mean(normAx(ii,:));'])
                eval([phaseNames{j},'.AyMean(ii) = mean(normAy(ii,:));'])            
                eval([phaseNames{j},'.PxSE(ii) = std(normPx(ii,:))/sqrt(length(normPx(ii,:)));'])
                eval([phaseNames{j},'.PySE(ii) = std(normPy(ii,:))/sqrt(length(normPy(ii,:)));'])
                eval([phaseNames{j},'.VxSE(ii) = std(normVx(ii,:))/sqrt(length(normVx(ii,:)));'])
                eval([phaseNames{j},'.VySE(ii) = std(normVy(ii,:))/sqrt(length(normVy(ii,:)));'])
                eval([phaseNames{j},'.FxSE(ii) = std(normFx(ii,:))/sqrt(length(normFx(ii,:)));'])
                eval([phaseNames{j},'.FySE(ii) = std(normFy(ii,:))/sqrt(length(normFy(ii,:)));'])
                eval([phaseNames{j},'.AxSE(ii) = std(normAx(ii,:))/sqrt(length(normAx(ii,:)));'])
                eval([phaseNames{j},'.AySE(ii) = std(normAy(ii,:))/sqrt(length(normAy(ii,:)));'])
            end
            %%% SUBJECT DATA
            eval(['numSubjs = size(',phaseNames{j},'.PxSubjNorm,2);'])
    
            %FOR NORMED SUBJECT MEAN 
            eval(['trLengthSubj = ones(1,numSubjs)*size(',phaseNames{j},'.PxSubjNorm,1);'])
            for jj = 1:numSubjs
                eval(['notLongest = ~isempty(find(isnan(',phaseNames{j},'.PxSubjNorm(:,jj)),1,''first''));'])
                if notLongest    
                    eval(['trLengthSubj(jj) = find(isnan(',phaseNames{j},'.PxSubjNorm(:,jj)),1,''first'')-1;'])
                end
            end
            avgTrLengthSubj = floor(mean(trLengthSubj));
        
        
            % Resample Subject Trajectories so the same time length
            normSubjPx = zeros(avgTrLengthSubj,numSubjs);
            normSubjPy = zeros(avgTrLengthSubj,numSubjs);
            normSubjVx = zeros(avgTrLengthSubj,numSubjs);
            normSubjVy = zeros(avgTrLengthSubj,numSubjs);
            normSubjFx = zeros(avgTrLengthSubj,numSubjs);
            normSubjFy = zeros(avgTrLengthSubj,numSubjs);
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
                eval(['normSubjFx(:,jj) = resample(',phaseNames{j},...
                    '.FxSubjNorm(1:trLengthSubj(jj),jj),avgTrLengthSubj,trLengthSubj(jj),0);'])            
                eval(['normSubjFy(:,jj) = resample(',phaseNames{j},...
                    '.FySubjNorm(1:trLengthSubj(jj),jj),avgTrLengthSubj,trLengthSubj(jj),0);'])
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
            eval([phaseNames{j},'.normSubjFx = normSubjFx;'])
            eval([phaseNames{j},'.normSubjFy = normSubjFy;'])
            eval([phaseNames{j},'.normSubjAx = normSubjAx;'])
            eval([phaseNames{j},'.normSubjAy = normSubjAy;'])
        
            % Find the average trajectory
            for ii = 1:avgTrLengthSubj 
                eval([phaseNames{j},'.PxSubjMean(ii) = mean(normSubjPx(ii,:));'])
                eval([phaseNames{j},'.PySubjMean(ii) = mean(normSubjPy(ii,:));'])
                eval([phaseNames{j},'.VxSubjMean(ii) = mean(normSubjVx(ii,:));'])
                eval([phaseNames{j},'.VySubjMean(ii) = mean(normSubjVy(ii,:));'])
                eval([phaseNames{j},'.FxSubjMean(ii) = mean(normSubjFx(ii,:));'])
                eval([phaseNames{j},'.FySubjMean(ii) = mean(normSubjFy(ii,:));'])    
                eval([phaseNames{j},'.AxSubjMean(ii) = mean(normSubjAx(ii,:));'])
                eval([phaseNames{j},'.AySubjMean(ii) = mean(normSubjAy(ii,:));'])            
                eval([phaseNames{j},'.PxSubjSE(ii) = std(normSubjPx(ii,:))/sqrt(length(normSubjPx(ii,:)));'])
                eval([phaseNames{j},'.PySubjSE(ii) = std(normSubjPy(ii,:))/sqrt(length(normSubjPy(ii,:)));'])
                eval([phaseNames{j},'.VxSubjSE(ii) = std(normSubjVx(ii,:))/sqrt(length(normSubjVx(ii,:)));'])
                eval([phaseNames{j},'.VySubjSE(ii) = std(normSubjVy(ii,:))/sqrt(length(normSubjVy(ii,:)));'])
                eval([phaseNames{j},'.FxSubjSE(ii) = std(normSubjFx(ii,:))/sqrt(length(normSubjFx(ii,:)));'])
                eval([phaseNames{j},'.FySubjSE(ii) = std(normSubjFy(ii,:))/sqrt(length(normSubjFy(ii,:)));'])
                eval([phaseNames{j},'.AxSubjSE(ii) = std(normSubjAx(ii,:))/sqrt(length(normSubjAx(ii,:)));'])
                eval([phaseNames{j},'.AySubjSE(ii) = std(normSubjAy(ii,:))/sqrt(length(normSubjAy(ii,:)));'])
            end
                
            % Calculate Standard Error
            eval([phaseNames{j},'.estGainMean = nanmean(',phaseNames{j},'.estGain);'])
            eval([phaseNames{j},'.estGainSE = nanstd(',phaseNames{j},'.estGain)/sqrt(sum(~isnan(',phaseNames{j},'.estGain)));'])
    
            % Max Horizontal Force
            eval([phaseNames{j},'.FxMaxMean = nanmean(',phaseNames{j},'.FxMax);'])
            eval([phaseNames{j},'.FxMaxSE = nanstd(',phaseNames{j},'.FxMax)/sqrt(sum(~isnan(',phaseNames{j},'.FxMax)));'])
             
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
            eval([phaseNames{j},'.normFx = [];'])
            eval([phaseNames{j},'.normFy = [];'])  
            eval([phaseNames{j},'.normAx = [];'])
            eval([phaseNames{j},'.normAy = [];'])            
            eval([phaseNames{j},'.PxMean = [];'])
            eval([phaseNames{j},'.PyMean = [];'])
            eval([phaseNames{j},'.VxMean = [];'])
            eval([phaseNames{j},'.VyMean = [];'])
            eval([phaseNames{j},'.FxMean = [];'])
            eval([phaseNames{j},'.FyMean = [];'])    
            eval([phaseNames{j},'.AxMean = [];'])
            eval([phaseNames{j},'.AyMean = [];'])            
            eval([phaseNames{j},'.PxSE = [];'])
            eval([phaseNames{j},'.PySE = [];'])
            eval([phaseNames{j},'.VxSE = [];'])
            eval([phaseNames{j},'.VySE = [];'])
            eval([phaseNames{j},'.FxSE = [];'])
            eval([phaseNames{j},'.FySE = [];'])
            eval([phaseNames{j},'.AxSE = [];'])
            eval([phaseNames{j},'.AySE = [];'])
            %
            eval([phaseNames{j},'.normSubjPx = [];'])
            eval([phaseNames{j},'.normSubjPy = [];'])
            eval([phaseNames{j},'.normSubjVx = [];'])
            eval([phaseNames{j},'.normSubjVy = [];'])
            eval([phaseNames{j},'.normSubjFx = [];'])
            eval([phaseNames{j},'.normSubjFy = [];'])   
            eval([phaseNames{j},'.normSubjAx = [];'])
            eval([phaseNames{j},'.normSubjAy = [];'])            
            eval([phaseNames{j},'.PxSubjMean = [];'])
            eval([phaseNames{j},'.PySubjMean = [];'])
            eval([phaseNames{j},'.VxSubjMean = [];'])
            eval([phaseNames{j},'.VySubjMean = [];'])
            eval([phaseNames{j},'.FxSubjMean = [];'])
            eval([phaseNames{j},'.FySubjMean = [];'])     
            eval([phaseNames{j},'.AxSubjMean = [];'])
            eval([phaseNames{j},'.AySubjMean = [];'])            
            eval([phaseNames{j},'.PxSubjSE = [];'])
            eval([phaseNames{j},'.PySubjSE = [];'])
            eval([phaseNames{j},'.VxSubjSE = [];'])
            eval([phaseNames{j},'.VySubjSE = [];'])
            eval([phaseNames{j},'.FxSubjSE = [];'])
            eval([phaseNames{j},'.FySubjSE = [];'])
            eval([phaseNames{j},'.AxSubjSE = [];'])
            eval([phaseNames{j},'.AySubjSE = [];'])
            %
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
    eval([phaseNames{ii},'.FxMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.FyMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AxMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AyMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.PxSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.PySENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.VxSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.VySENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.FxSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.FySENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AxSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AySENonNorm = zeros(1,maxTrialLength);']);
    for jj = 1:maxTrialLength
        eval([phaseNames{ii},'.PxMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Px(jj,:));'])
        eval([phaseNames{ii},'.PyMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Py(jj,:));'])
        eval([phaseNames{ii},'.VxMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Vx(jj,:));'])
        eval([phaseNames{ii},'.VyMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Vy(jj,:));'])
        eval([phaseNames{ii},'.FxMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Fx(jj,:));'])
        eval([phaseNames{ii},'.FyMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Fy(jj,:));'])
        eval([phaseNames{ii},'.AxMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Ax(jj,:));'])
        eval([phaseNames{ii},'.AyMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.Ay(jj,:));'])
        
        % Commented Out to Inspect STD instead of SE
        eval([phaseNames{ii},'.PxSENonNorm(jj) = nanstd(',phaseNames{ii},'.Px(jj,:))/sqrt(nTrials);'])
        eval([phaseNames{ii},'.PySENonNorm(jj) = nanstd(',phaseNames{ii},'.Py(jj,:))/sqrt(nTrials);'])
        eval([phaseNames{ii},'.VxSENonNorm(jj) = nanstd(',phaseNames{ii},'.Vx(jj,:))/sqrt(nTrials);'])
        eval([phaseNames{ii},'.VySENonNorm(jj) = nanstd(',phaseNames{ii},'.Vy(jj,:))/sqrt(nTrials);'])
        eval([phaseNames{ii},'.FxSENonNorm(jj) = nanstd(',phaseNames{ii},'.Fx(jj,:))/sqrt(nTrials);'])
        eval([phaseNames{ii},'.FySENonNorm(jj) = nanstd(',phaseNames{ii},'.Fy(jj,:))/sqrt(nTrials);'])
        eval([phaseNames{ii},'.AxSENonNorm(jj) = nanstd(',phaseNames{ii},'.Ax(jj,:))/sqrt(nTrials);'])
        eval([phaseNames{ii},'.AySENonNorm(jj) = nanstd(',phaseNames{ii},'.Ay(jj,:))/sqrt(nTrials);'])

    end
    
    % Truncate to mean length.
    eval(['meanLength =  length(',phaseNames{ii},'.PxMean);']);
    eval([phaseNames{ii},'.PxMeanNonNorm = ',phaseNames{ii},'.PxMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.PyMeanNonNorm = ',phaseNames{ii},'.PyMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VxMeanNonNorm = ',phaseNames{ii},'.VxMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VyMeanNonNorm = ',phaseNames{ii},'.VyMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.FxMeanNonNorm = ',phaseNames{ii},'.FxMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.FyMeanNonNorm = ',phaseNames{ii},'.FyMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.AxMeanNonNorm = ',phaseNames{ii},'.AxMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.AyMeanNonNorm = ',phaseNames{ii},'.AyMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.PxSENonNorm = ',phaseNames{ii},'.PxSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.PySENonNorm = ',phaseNames{ii},'.PySENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VxSENonNorm = ',phaseNames{ii},'.VxSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VySENonNorm = ',phaseNames{ii},'.VySENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.FxSENonNorm = ',phaseNames{ii},'.FxSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.FySENonNorm = ',phaseNames{ii},'.FySENonNorm(1:meanLength);']);
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
    eval([phaseNames{ii},'.FxSubjMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.FySubjMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AxSubjMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AySubjMeanNonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.PxSubjSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.PySubjSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.VxSubjSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.VySubjSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.FxSubjSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.FySubjSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AxSubjSENonNorm = zeros(1,maxTrialLength);']);
    eval([phaseNames{ii},'.AySubjSENonNorm = zeros(1,maxTrialLength);']);
    for jj = 1:maxTrialLength
        eval([phaseNames{ii},'.PxSubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.PxSubj(jj,:));'])
        eval([phaseNames{ii},'.PySubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.PySubj(jj,:));'])
        eval([phaseNames{ii},'.VxSubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.VxSubj(jj,:));'])
        eval([phaseNames{ii},'.VySubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.VySubj(jj,:));'])
        eval([phaseNames{ii},'.FxSubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.FxSubj(jj,:));'])
        eval([phaseNames{ii},'.FySubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.FySubj(jj,:));'])
        eval([phaseNames{ii},'.AxSubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.AxSubj(jj,:));'])
        eval([phaseNames{ii},'.AySubjMeanNonNorm(jj) = nanmean(',phaseNames{ii},'.AySubj(jj,:));'])
        
        % Commented Out to Inspect STD instead of SE
        eval([phaseNames{ii},'.PxSubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.PxSubj(jj,:))/sqrt(nSubj);'])
        eval([phaseNames{ii},'.PySubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.PySubj(jj,:))/sqrt(nSubj);'])
        eval([phaseNames{ii},'.VxSubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.VxSubj(jj,:))/sqrt(nSubj);'])
        eval([phaseNames{ii},'.VySubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.VySubj(jj,:))/sqrt(nSubj);'])
        eval([phaseNames{ii},'.FxSubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.FxSubj(jj,:))/sqrt(nSubj);'])
        eval([phaseNames{ii},'.FySubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.FySubj(jj,:))/sqrt(nSubj);'])
        eval([phaseNames{ii},'.AxSubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.AxSubj(jj,:))/sqrt(nSubj);'])
        eval([phaseNames{ii},'.AySubjSENonNorm(jj) = nanstd(',phaseNames{ii},'.AySubj(jj,:))/sqrt(nSubj);'])
    end
      
    if age == 'Y' || age == 'y'
        trialLengths = [143,114,96,111,111,110,99,102];
    else
        trialLengths = [148,172,126,120,139,129,119,133]; %TBD
    end
    meanLength = trialLengths(ii);
    
    eval([phaseNames{ii},'.PxSubjMeanNonNorm = ',phaseNames{ii},'.PxSubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.PySubjMeanNonNorm = ',phaseNames{ii},'.PySubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VxSubjMeanNonNorm = ',phaseNames{ii},'.VxSubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VySubjMeanNonNorm = ',phaseNames{ii},'.VySubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.AxSubjMeanNonNorm = ',phaseNames{ii},'.AxSubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.AySubjMeanNonNorm = ',phaseNames{ii},'.AySubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.FxSubjMeanNonNorm = ',phaseNames{ii},'.FxSubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.FySubjMeanNonNorm = ',phaseNames{ii},'.FySubjMeanNonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.PxSubjSENonNorm = ',phaseNames{ii},'.PxSubjSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.PySubjSENonNorm = ',phaseNames{ii},'.PySubjSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VxSubjSENonNorm = ',phaseNames{ii},'.VxSubjSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.VySubjSENonNorm = ',phaseNames{ii},'.VySubjSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.FxSubjSENonNorm = ',phaseNames{ii},'.FxSubjSENonNorm(1:meanLength);']);
    eval([phaseNames{ii},'.FySubjSENonNorm = ',phaseNames{ii},'.FySubjSENonNorm(1:meanLength);']);
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