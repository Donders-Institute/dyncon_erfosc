function concentric_grating_experiment(fileName, isLive, drawMask)
% This function is copied from DriftDemo and modified to display a (masked)
% animated concentric grating moving inward. At random time points, the
% grating will make a 180 degree phase shift, which has to be reported by
% the subject by a button press. 20% of the trials have no phase shift.

% NOTE: costs a lot of memory because all images are made individually

% _________________________________________________________________________
%
% This is a very simple, bare bones demo on how to do frame animation. For
% much more efficient ways to draw gratings and gabors, have a look at
% DriftDemo2, DriftDemo3, DriftDemo4, ProceduralGaborDemo, GarboriumDemo,
% ProceduralGarboriumDemo and DriftWaitDemo.

% NOTE: DriftDemo2 is more efficient because only one texture is
% generated, which is moved every frame to create the illusion of movement.
% However, this is not possible for a concentric grating, so every frame a
% new texture has to be generated.
% _________________________________________________________________________

% HISTORY
% 7/6/2016: change grating to concentric grating
% 8/6/2016: add circular and gaussian mask, introduce phase-shift
% 9/6/2016: make trials, introduce baseline period
% 10/6/2016: add no shift trials
% 15/5/2016: change background color
% 16/6/2016: edit Gaussian mask, add (rectangular) hanning mask
% 16/6/2016: add circular hanning mask

%% Function settings


if nargin<1
    fileName = 'test';
end

if isempty(fileName)
    fileName = 'test';
end

if nargin<2
    isLive = 0;
end

if isempty(isLive)
    isLive = 0;
end

if nargin<3
    drawMask = 1;
end

if isempty(drawMask)
    drawMask = 1;
end
try
    Screen('Preference', 'SkipSyncTests', 1); % THIS REMOVES THE SYNCHRONIZATION CHECK! REMOVE BEFORE EXPERIMENT
    %% Standard PTB settings
    % This script calls Psychtoolbox commands available only in OpenGL-based
    % versions of the Psychtoolbox. (So far, the OS X Psychtoolbox is the
    % only OpenGL-base Psychtoolbox.)  The Psychtoolbox command AssertPsychOpenGL will issue
    % an error message if someone tries to execute this script on a computer without
    % an OpenGL Psychtoolbox
    AssertOpenGL;
    
    % Get the list of screens and choose the one with the highest screen number.
    % Screen 0 is, by definition, the display with the menu bar. Often when
    % two monitors are connected the one without the menu bar is used as
    % the stimulus display.  Chosing the display with the highest dislay number is
    % a best guess about where you want the stimulus displayed.
    screens=Screen('Screens');
    if isLive
        screenNumber=max(screens);
    else
        screenNumber=2;
    end
    
    % Find the color values which correspond to white and black: Usually
    % black is always 0 and white 255, but this rule is not true if one of
    % the high precision framebuffer modes is enabled via the
    % PsychImaging() commmand, so we query the true values via the
    % functions WhiteIndex and BlackIndex:
    white=WhiteIndex(screenNumber);
    black=BlackIndex(screenNumber);
    
    % Round gray to integral number, to avoid roundoff artifacts with some
    % graphics cards:
    gray=round((white+black)/2);
    
    % This makes sure that on floating point framebuffers we still get a
    % well defined gray. It isn't strictly neccessary in this demo:
    if gray == white
        gray=white / 2;
    end
    
    % Contrast 'inc'rement range for given white and gray values:
    inc=white-gray;
    
    % Open a double buffered fullscreen window and select a black background
    % color:
    [window, windowRect]=Screen('OpenWindow',screenNumber, gray);
    
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(windowRect);
    
    % Query the frame duration
    ifi = Screen('GetFlipInterval', window);
    frameRate=Screen('FrameRate',screenNumber); %1/ifi
    
    % Set up alpha-blending for smooth (anti-aliased) lines
    Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
    
    %     HideCursor;
    
    
    %% Trigger settings
    if isLive
        btsi = Bitsi('com1');% when there is no com connected, use keyboard as input
    else
        btsi = Bitsi('');
    end
    
    xp = expconstants();
    
    %% Input settings
    
    KbName('UnifyKeyNames'); % converts operating specific keynames to universal keyname
    escape = KbName('ESCAPE');
    space = KbName('space');
    
    %% Experimental settings
    % calculate the size of the stimulus according to the used setup
    if isLive
        screen.width = 49; %MEG screenwidth in cm
        screen.viewingDistance = 76; % in cm
        screen.resolutionX = 1024; %MEG screen in pixels
    else
        screen.width = 52;%desktop screenwidth in cm
        screen.viewingDistance = 59; % in cm
        screen.resolutionX = 1920; %desktop screen in pixels
    end
    
    screen.totVisDeg = 2*atan(screen.width/(2*screen.viewingDistance))*(180/pi); % formula for calculating visual degrees
    screen.pixPerDeg = screen.resolutionX/screen.totVisDeg;
    
    visualAngleGrating = 12.5;
    gratingSize = visualAngleGrating*screen.pixPerDeg; % calculate how big the grating should be (can not be rounded off, might result in odd integer)
    gratingRadius = round(gratingSize/2); % the grating can only exist of integers, so round
    gratingSize = 2*gratingRadius; % to prevent consistency errors, redifine gratingSize
    
    instructions = 'Instructions \n\n Fixate at the middle of the screen. \n\n Press the button if you see a shift in the animation';
    ntrials=200;
    chanceNoShift = 0.2; % in 20% of the trials, have no shift.
    preStimTime=2; % fixation, baseline period
    preShiftTime=1.5; % cue can come 1.5 seconds after start trial
    shiftRange = 5; % present the shift at random in 5 seconds
    postShiftTime = 1.5; % keep animation going 1.5 seconds after shift
    driftFreq = 1; % in one second, every pixel of the grating completes one cylce (black-white-black)
    numFrames=round(driftFreq/ifi); % temporal period, in frames, of the drifting grating
    
    % settings grating
    [x,y]=meshgrid(-gratingRadius:gratingRadius,-gratingRadius:gratingRadius);
    f=1.5*2*pi; % cycles/pixel
    
    keyPress = zeros(ntrials,4); %safe shift onset, response and response time
    
    %% Generate stimulus
    % generate a mask to cover the edges of the grating, making it a circle
    %     R = gratingRadius ; % radius
    %     M = zeros(2*gratingRadius+1, 2*gratingRadius+1, 1); %zero matrix same size as the grid
    %     M(:,:,1) = sqrt(x.^2 + y.^2) < R ; %disk formula: all masked parts of the grid must be zero.
    
    % make hanning mask
    L = 2*gratingRadius+1;
    w1D = hann(L); % 1D hann window
    xx = linspace(-gratingRadius,gratingRadius,L);
    [X,Y] = meshgrid(xx);
    r = sqrt( X.^2 + Y.^2 );
    w2D = zeros(L);
    % the hanning window has to go from 1 to 0.5
    w2D(r<=gratingRadius) = 0.5*interp1(xx,w1D,r(r<=gratingRadius)); % first go from 0.5 to 0
    w2D = w2D + 0.5; %add 0.5 to go from 1 to 0.5
    
    %%
    % Generate grating texture
    % Compute each frame of the movie and convert those frames stored in
    % MATLAB matrices, into Psychtoolbox OpenGL textures using 'MakeTexture';
    % NOTE: very inefficient, costs a lot of memory.. alternative?
    tex=zeros(numFrames,1);
    
    for i=1:numFrames
        phase=(i/numFrames)*2*pi; %change the phase of the grating according to the framenumber
        m=sin(sqrt(x.^2+y.^2) / f + phase);
        grating = w2D.*(gray+inc*m); %multiply the hanning mask with the
        % grating. gray takes care of coloring. inc*m is the grating
        % itself.
        grating(sqrt(x.^2+y.^2)>gratingRadius) = gray; % make everything
        % outside the circle the same color as the background (gray)
        tex(i)=Screen('MakeTexture', window, grating); % multiply mask with grating
    end
       
    
    % _____________________________________________________________________
    % Other masks
    %{

    % w2D(:, :, 2)=white * (1 - 0.5*exp(-((x/120).^2)-((y/120).^2))); % Gaussian going down to gray (for black remove 0.5*)
%     for n=1:2*gratingRadius+1
%     w2D(:,n,2) = hann(2* gratingRadius+1); % apply 1D hanning over the columns
%     end
%     w2D(:,:,2) = gray*mask(:,:,1).* ( 1 - (mask(:,:,2).*mask(:,:,2)')); % multiply with
    % inverse to apply hanning window to rows as well. Invert (1-hanning)
    % to get lowest value in centre. Multiply with white to get the right colorscale.
    % NOTE: This Hanning window is rectangular, not circular! That is why
    % the window is multiplied with a circular mask (mask(:,:,1))
      masktex=Screen('MakeTexture', window, w2D);

    %}
    % _____________________________________________________________________
    
    shiftPoint = ((preShiftTime+shiftRange)-preShiftTime).*rand(ntrials,1) + preShiftTime;
    trlNoShift = sort(randperm(ntrials, ntrials*chanceNoShift)); % generate a number of random (w/o replacement) trials that won't have a phase shift
    %   shiftPoint(trlNoShift)=0; % replace the noShift trials with 0.
    %   instead, make sure that these trials won't have a shift (otherwise it
    %   also means short trial)
    shiftFrame = round(shiftPoint*frameRate);
    keyPress(:,1)=shiftPoint;
    
    KbName('UnifyKeyNames'); % convert operating specific keynames to general keynames
    escape = KbName('ESCAPE');
    space = KbName('space');
    
    %% Start Experiment
    
    % Run the movie animation for a variable period: end it 1.5 seconds
    % after phase shift.
    movieDurationSecs=shiftPoint+postShiftTime;
    % Convert movieDuration in seconds to duration in frames to draw:
    movieDurationFrames=round(movieDurationSecs * frameRate);
    % assign the right texture indices to frames
    for trl=1:ntrials
        if trl==1
            % Provide instructions at the beginning of the experiment
            Screen('FillRect', window, gray);
            DrawFormattedText(window, instructions,...
                'center', 'center', white, 55, [], [], 1.3);
            Screen('Flip', window);
            KbWait; % wait until key is pressed
            pause(0.5)
        end
        movieFrameIndices=mod(0:(movieDurationFrames(trl)-1), numFrames) + 1; %assign the right image index to each frame
        
        
        % Use realtime priority for better timing precision:
        priorityLevel=MaxPriority(window);
        Priority(priorityLevel);
        
        % Baseline period, black screen with fixation dot
        Screen('DrawDots', window, [xCenter yCenter], 20, [white white white], [], 2);
        Screen('Flip', window)
        pause(preStimTime)
        tic % measure time
        
        btsi.sendTrigger(xp.TRIG_ONSET_TRIAL);
        % Animation loop:
        for i=1:movieDurationFrames(trl)
            if i==shiftFrame(trl)
                t0=GetSecs; % get time stamp at the moment of phase shift in order to calculate reaction time
            end
            
            % at the moment the phase-shift should occur, jump to the texture
            % with a 180 degree phase shift, and add this shift for the
            % remainder of the video. Otherwise, don't shift.
            if i>=shiftFrame(trl) && movieFrameIndices(i)>numFrames/2 && ~any(trlNoShift==trl)
                shift = -numFrames/2; % if more than half of the period already passed, subtract the shift
                % enable abort using esc
            elseif i>=shiftFrame(trl) && movieFrameIndices(i)<=numFrames/2 && ~any(trlNoShift==trl)
                shift = numFrames/2; % otherwise add the shift (because the period loops)
            else
                shift=0;
            end
            
            % Draw the appropriate image
            Screen('DrawTexture', window, tex(movieFrameIndices(i) + shift));
            
            % Draw the Gaussian mask if required
            if drawMask
                %                 Screen('DrawTexture', win dow, masktex);
            end
            
            % Show it at next display vertical retrace. Please check DriftDemo2
            % and later, as well as DriftWaitDemo for much better approaches to
            % guarantee a robust and constant animation display timing! This is
            % very basic and not best practice!
            if i==shiftFrame(trl);
                btsi.sendTrigger(xp.TRIG_SHIFT);
            end
            Screen('Flip', window);
            [keyIsDown, secs, keyCode] = KbCheck(-3);
            %% How to make sure only the first response is saved?
            %             while keyCode(space)==0
            %                 [keyIsDown, secs, keyCode] = KbCheck(-3); % check whether a button was pressed; only record first button press
            %             end
            if keyCode(escape) %if escape is pressed
                save(fileName, 'keyPress', 'trlNoShift') % save the data
                sca % close the screen
                break % and end PTB
            elseif keyCode(space) % if space is pressed
                if i>=shiftFrame(trl)
                    keyPress(trl,2)=1; % report
                    keyPress(trl,3) = secs-t0; % report reaction time
                else
                    keyPress(trl,2)=2; % if space is pressed before phase shift, report it.
                end
            end
            
        end
        btsi.sendTrigger(xp.TRIG_OFFSET_TRIAL);
        keyPress(trl,4)=toc-postShiftTime;% the amount of time it takes for the animation
        % keyPress(:,1) is the amount of time it should take to run the
        % animation.
        Priority(0);
    end
    save(fileName, 'keyPress', 'trlNoShift')
    
    
    % Close all textures. This is not strictly needed, as
    % Screen('CloseAll') would do it anyway. However, it avoids warnings by
    % Psychtoolbox about unclosed textures. The warnings trigger if more
    % than 10 textures are open at invocation of Screen('CloseAll') and we
    % have 12 textues here:
    Screen('Close');
    
    % Close window:
    Screen('CloseAll');
    
catch
    %this "catch" section executes in case of an error in the "try" section
    %above.  Importantly, it closes the onscreen window if its open.
    Screen('CloseAll');
    Priority(0);
    psychrethrow(psychlasterror);
end %try..catch..
