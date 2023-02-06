%% Stats on Thrown Out Trials


% By Age Group
ages = {'Y','O'};

% By Subject Number
subj_num = [13,15];

% For Each Phase
% phaseNames = {'EN1','LN1','EF1','LF1','EF2','LF2','EN2','LN2'};
phaseNames = {'LN1','EF1','LF2','EN2'};

tots{length(ages)} = [];

for ii = 1:length(ages)
    if ages{ii} == 'Y'
        load('../../../Data/Trajectory Data/Processed Data/Younger-AveragedTrajectoryData.mat','yData','yDataChan')
    else
        load('../../../Data/Trajectory Data/Processed Data/Older-AveragedTrajectoryData.mat','yData','yDataChan')
    end

    tot_bad = 0;
    tot_fast = 0;
    tot_slow = 0;
    tot = 0;

    tot_bad_chan = 0;
    tot_fast_chan = 0;
    tot_slow_chan = 0;
    tot_chan = 0;

    tot_bad_all = [];
    tot_all = [];
    tot_bad_all_chan = [];
    tot_all_chan = [];
    
    
    for jj = 1:subj_num(ii)
        for kk = 1:length(phaseNames)
            chans = [];
            tot_bad_all(jj,kk) = 0;
            tot_all(jj,kk) = 0;
            tot_bad_all_chan(jj,kk) = 0;
            tot_all_chan(jj,kk) = 0;  
            
           eval(['tot_bad = tot_bad + length(yData.',phaseNames{kk},...
               '.trialNums{jj}.tooLong) + length(yData.',phaseNames{kk},...
               '.trialNums{jj}.tooShort);']);
           eval(['tot_fast = tot_fast + length(yData.',phaseNames{kk},...
               '.trialNums{jj}.tooShort);']);
           eval(['tot_slow = tot_slow + length(yData.',phaseNames{kk},...
               '.trialNums{jj}.tooLong);']);
           eval(['tot = tot + length(yData.',phaseNames{kk},...
               '.trialNums{jj}.tooLong) + length(yData.',phaseNames{kk},...
               '.trialNums{jj}.tooShort) + length(yData.',phaseNames{kk},...
               '.trialNums{jj}.justRight);']);
            % Channel
           eval(['tot_bad_chan = tot_bad_chan + length(yDataChan.',phaseNames{kk},...
               '.ctrialNums{jj}.tooLong) + length(yDataChan.',phaseNames{kk},...
               '.ctrialNums{jj}.tooShort);']);
           eval(['tot_fast_chan = tot_fast_chan + length(yDataChan.',phaseNames{kk},...
               '.ctrialNums{jj}.tooShort);']);
           eval(['tot_slow_chan = tot_slow_chan + length(yDataChan.',phaseNames{kk},...
               '.ctrialNums{jj}.tooLong);']);
           eval(['tot_chan = tot_chan + length(yDataChan.',phaseNames{kk},...
               '.ctrialNums{jj}.tooLong) + length(yDataChan.',phaseNames{kk},...
               '.ctrialNums{jj}.tooShort) + length(yDataChan.',phaseNames{kk},...
               '.ctrialNums{jj}.justRight);']);

           % Totally Separated
           eval(['tot_bad_all(jj,kk) = tot_bad_all(jj,kk) + length(yData.',phaseNames{kk},...
               '.trialNums{jj}.tooLong) + length(yData.',phaseNames{kk},...
               '.trialNums{jj}.tooShort);']);
           eval(['tot_all(jj,kk) = tot_all(jj,kk) + length(yData.',phaseNames{kk},...
               '.trialNums{jj}.tooLong) + length(yData.',phaseNames{kk},...
               '.trialNums{jj}.tooShort) + length(yData.',phaseNames{kk},...
               '.trialNums{jj}.justRight);']);        
           eval(['tot_bad_all_chan(jj,kk) = tot_bad_all_chan(jj,kk) + length(yDataChan.',phaseNames{kk},...
               '.ctrialNums{jj}.tooLong) + length(yDataChan.',phaseNames{kk},...
               '.ctrialNums{jj}.tooShort);']);
           eval(['tot_all_chan(jj,kk) = tot_all_chan(jj,kk) + length(yDataChan.',phaseNames{kk},...
               '.ctrialNums{jj}.tooLong) + length(yDataChan.',phaseNames{kk},...
               '.ctrialNums{jj}.tooShort) + length(yDataChan.',phaseNames{kk},...
               '.ctrialNums{jj}.justRight);']);
        end
    end
    tots{ii}.tot_bad_all = tot_bad_all;
    tots{ii}.tot_all = tot_all;
    tots{ii}.tot_bad_all_chan = tot_bad_all_chan;
    tots{ii}.tot_all_chan = tot_all_chan;
    tot_bad
    tot_fast
    tot_slow
    tot
    tot_bad_chan
    tot_fast_chan
    tot_slow_chan
    tot_chan
end


for ii = 1:length(ages)
    if ages{ii} == 'Y'
        load('../../../Data/Trajectory Data/Processed Data/Younger-AveragedTrajectoryData.mat','yData','yDataChan')
    else
        load('../../../Data/Trajectory Data/Processed Data/Older-AveragedTrajectoryData.mat','yData','yDataChan')
    end

    for kk = 1:length(phaseNames)
        temp = [];
        for jj = 1:subj_num(ii)
%             eval(['disp(num2str(yDataChan.',phaseNames{kk},'.ctrialNums{jj}.justRight))'])
            eval(['exist =  ~isempty(yDataChan.',phaseNames{kk},'.ctrialNums{jj}.justRight);'])
            if exist
                eval(['temp = [temp, yDataChan.',phaseNames{kk},'.ctrialNums{jj}.justRight''];'])
            end
        end
        chans{ii,kk} = temp;
    end
end