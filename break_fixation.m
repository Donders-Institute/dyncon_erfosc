function [sample_breakFixation] = break_fixation(data, visDegOffFixation)

for iTrial = 1:size(data.trial,2)
    isFixation = sqrt(visAngleX{iTrial}.^2 + visAngleY{iTrial}.^2) < visDegOffFixation; % note all the samples where fixation is good/bad
    %         samplenumber = data.sampleinfo(iTrial,1);
    %         while samplenumber <= data.sampleinfo(iTrial,2)
    samplenumber = 1;
    while samplenumber <= size(isFixation,2) % for the whole trial (in sample numbers)
        % if there is a new off-fixation event, note the first off-fixation sample
        if any(isFixation(1, samplenumber:end)==0)
            sampleStartOffFixation = (samplenumber -1) + find(isFixation(1, samplenumber:end)==0,1);
            sample_breakFixation(end+1,1) = data.sampleinfo(iTrial,1) + sampleStartOffFixation;
            %                 samplenumber = sampleStartOffFixation+1;
            %if there is fixation returns after off-fixation, note last
            % off-fixation sample of that event
            if any(isFixation(1, sampleStartOffFixation + 1 : end)==1)
                sampleStopOffFixation = sampleStartOffFixation + find(isFixation(1, sampleStartOffFixation + 1 : end) ==1, 1) -1; % stop off-fixation is 1 sample before on fixation again.
                sample_breakFixation(end,2) = data.sampleinfo(iTrial,1) + sampleStopOffFixation;
                samplenumber = sampleStopOffFixation +1;
break
            else
                sample_breakFixation(end,2) = data.sampleinfo(iTrial,2);
                samplenumber = data.sampleinfo(iTrial,2);
            end
        else
            break
        end
    end
end